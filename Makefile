EXTRA_CXXFLAGS=
CXXFLAGS=-ggdb -O3 -Wall -std=c++17 $(EXTRA_CXXFLAGS)

all: uvid_compress uvid_decompress

uvid_compress: uvid_compress.o
	g++ -o $@ $^

uvid_decompress: uvid_decompress.o
	g++ -o $@ $^



clean:
	rm -f uvid_compress uvid_decompress *.o
