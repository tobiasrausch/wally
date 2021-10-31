# use the ubuntu base image
FROM ubuntu:18.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install wally
RUN cd /opt \
    && git clone --recursive https://github.com/tobiasrausch/wally.git \
    && cd /opt/wally/ \
    && make STATIC=1 all \
    && make install


# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/wally/bin
WORKDIR /opt/wally/bin
COPY --from=0 /opt/wally/bin/wally .

# Workdir
WORKDIR /root/

# Add Wally to PATH
ENV PATH="/opt/wally/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]