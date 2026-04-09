CONGA_VERSION := "1.0"
CONGA_UPDATE := "April 9, 2026"
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -DCONGA_VERSION=\"$(CONGA_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DCONGA_UPDATE=\"$(CONGA_UPDATE)\"
LDFLAGS = htslib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -no-pie
NOCRAMFLAGS = htslib/libhts.a -lz -lm -lpthread
SOURCES = svdepth.c cmdline.c common.c bam_data.c read_distribution.c free.c likelihood.c svs.c split_read.c genome.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = conga
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	make clean -C htslib
	rm -f $(EXECUTABLE) *.o *~

nocram: $(OBJECTS)
	cd htslib && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make && cd ..
	$(CC) $(OBJECTS) -o $(EXECUTABLE)-nocram $(NOCRAMFLAGS)

libs:
	make -C htslib

install:
	cp conga $(INSTALLPATH)
