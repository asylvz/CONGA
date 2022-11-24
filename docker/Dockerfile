FROM ubuntu:18.04

RUN apt-get update -y --fix-missing
RUN apt-get upgrade -y

RUN apt-get install git make gcc zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev -y

RUN mkdir /conga
RUN mkdir /input
RUN mkdir /output
WORKDIR /conga

RUN git clone --recursive https://github.com/asylvz/CONGA /conga
RUN make libs && make


WORKDIR /conga

RUN apt-get remove git -y
RUN apt-get autoremove -y
VOLUME /input
VOLUME /output
ENTRYPOINT ["/conga/conga"]
