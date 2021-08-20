# Dockerfile for Python based RDKit implementation
# Based on Debian.
# Includes InCHI support.
# WARNING: this takes about an hour to build

FROM debian:stretch
MAINTAINER Tim Dudgeon <tdudgeon@informaticsmatters.com>

RUN apt-get update && apt-get install -y \
 build-essential\
 python-numpy\
 cmake\
 python-dev\
 python-pip\
 sqlite3\
 libsqlite3-dev\
 libboost-dev\
 libboost-system-dev\
 libboost-thread-dev\
 libboost-serialization-dev\
 libboost-python-dev\
 libboost-regex-dev\
 swig\
 git\
 wget\
 zip &&\
 apt-get upgrade -y &&\
 apt-get clean -y
 
ENV RDKIT_BRANCH=master
RUN git clone -b $RDKIT_BRANCH --single-branch https://github.com/rdkit/rdkit.git

ENV RDBASE=/rdkit
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib:/usr/lib/x86_64-linux-gnu
ENV PYTHONPATH=$PYTHONPATH:$RDBASE

RUN mkdir $RDBASE/build
WORKDIR $RDBASE/build

RUN cmake -DRDK_BUILD_INCHI_SUPPORT=ON .. &&\
 make &&\
 make install &&\
 make clean

WORKDIR $RDBASE
uwsgi --http-socket 0.0.0.0:8080 --wsgi-file main.py --callable app --uid ubuntu
