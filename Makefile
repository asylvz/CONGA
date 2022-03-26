CONGA_VERSION := "0.1"
CONGA_UPDATE := "March 26, 2022"
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I sonic -DCONGA_VERSION=\"$(CONGA_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DCONGA_UPDATE=\"$(CONGA_UPDATE)\"
LDFLAGS = htslib/libhts.a sonic/libsonic.a -lz -lm -lpthread -llzma -lbz2 -lcurl -no-pie
SOURCES = svdepth.c cmdline.c common.c bam_data.c read_distribution.c free.c likelihood.c svs.c split_read.c
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
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib
	make -C sonic

install:
	cp Conga $(INSTALLPATH)
