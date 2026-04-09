CONGA_VERSION := "1.1"
CONGA_UPDATE := "April 9, 2026"
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I src -DCONGA_VERSION=\"$(CONGA_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DCONGA_UPDATE=\"$(CONGA_UPDATE)\"
LDFLAGS = htslib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl
NOCRAMFLAGS = htslib/libhts.a -lz -lm -lpthread
SOURCES = src/svdepth.c src/cmdline.c src/common.c src/bam_data.c src/read_distribution.c src/free.c src/likelihood.c src/svs.c src/split_read.c src/genome.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = conga
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf src/*.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	make clean -C htslib
	rm -f $(EXECUTABLE) src/*.o *~

nocram: $(OBJECTS)
	cd htslib && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make && cd ..
	$(CC) $(OBJECTS) -o $(EXECUTABLE)-nocram $(NOCRAMFLAGS)

libs:
	make -C htslib

install:
	cp conga $(INSTALLPATH)
