FROM ubuntu:22.04

RUN apt-get update -y --fix-missing && \
    apt-get install -y git make gcc zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev && \
    apt-get clean

RUN mkdir /input /output
WORKDIR /conga

RUN git clone --recursive https://github.com/asylvz/CONGA /conga && \
    make libs && make && \
    apt-get remove -y git && \
    apt-get autoremove -y

VOLUME /input
VOLUME /output
ENTRYPOINT ["/conga/conga"]
