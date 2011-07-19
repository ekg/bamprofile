CXX=g++
BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib
CXXFLAGS=-O3 -I$(BAMTOOLS_ROOT)/include -L./ #-L$(BAMTOOLS_ROOT)/lib

all: bamprofile

# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && make
	cp bamtools/lib/libbamtools.a ./

# statically compiles bamaddrg against bamtools static lib
bamprofile: bamprofile.cpp libbamtools.a
	$(CXX) $(CXXFLAGS) fastahack/Fasta.cpp fastahack/split.cpp bamprofile.cpp -o bamprofile -lbamtools -lz

clean:
	rm -f bamprofile libbamtools.a *.o

.PHONY: clean
