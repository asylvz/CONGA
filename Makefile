SVDEPTH_VERSION := "0.1"
SVDEPTH_UPDATE := "October 11, 2019"
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I sonic -DSVDEPTH_VERSION=\"$(SVDEPTH_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DSVDEPTH_UPDATE=\"$(SVDEPTH_UPDATE)\"
LDFLAGS = htslib/libhts.a sonic/libsonic.a -lz -lm -lpthread -llzma -lbz2 -lcurl
SOURCES = svdepth.c cmdline.c common.c bam_data.c read_distribution.c free.c likelihood.c svs.c split_read.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = svdepth
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib
	make -C sonic

install:
	cp SvDepth $(INSTALLPATH)
