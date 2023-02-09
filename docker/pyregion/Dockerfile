# syntax=docker/dockerfile:1
ARG TARGETPLATFORM=linux/amd64
FROM --platform=${TARGETPLATFORM} ubuntu:18.04

ARG TARGETPLATFORM
ENV TZ=Etc/UTC
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
    #   && apt-get install -yq tzdata \
    #   && ln -fs /usr/share/zoneinfo/${TZ} /etc/localtime \
    #   && dpkg-reconfigure -f noninteractive tzdata \
      && apt-get -y --no-install-recommends install \
         git \
         wget \
         python3-dev \
         python3-pip \
         python-pyparsing \
         build-essential \
        #  cmake \
        #  gdb \
        #  g++ \
        #  ca-certificates \
        #  gfortran \
        #  libgfortran3 \
        #  libc6 \
        #  libhwloc-dev \
        #  libibverbs-dev \
        #  libopenmpi-dev \
        #  imagemagick \
        #  libmagickcore-dev \
        #  libmagickwand-dev \
        #  libmagic-dev \
        #  fuse \
        #  libfuse2 \
        #  libclang-dev \
        #  libssl-dev \
        #  libffi-dev \
        #  libncurses5-dev \
        #  libreadline-dev \
        #  flex \
        #  bison \
        #  libblas-dev \
        #  liblapacke-dev \
        #  libcfitsio-dev \
        #  libpthread-stubs0-dev \
        #  libgsl-dev \
        #  wcslib-dev \
        #  libfftw3-dev \
         python-numpy \
         python-setuptools \
        #  libx11-dev \
        #  xauth \
        #  xorg \
         python3-pyqt5 \
         libboost-python-dev \
      && apt-get clean all \
      && rm -r /var/lib/apt/lists/*
WORKDIR /

# Update python3 to python 3.8 instead of the default python3.6
#RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y software-properties-common
#RUN add-apt-repository ppa:deadsnakes/ppa && DEBIAN_FRONTEND=noninteractive apt install -y python3.8
#RUN ln -s /usr/bin/pip3 /usr/bin/pip && \
#    ln -s /usr/bin/python3.8 /usr/bin/python

RUN : \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3.8-venv \
        python3.8-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && :

RUN python3.8 -m venv /venv
ENV PATH=/venv/bin:$PATH

# Set default python to python3:
#RUN update-alternatives --remove python /usr/bin/python2 \
#      && update-alternatives --install /usr/bin/python python /usr/bin/python3 10 \
#      && python3 -m pip install -U --force-reinstall pip \
#      && pip3 install setuptools \
#      && pip3 install wheel

# What version/branch of your repo do you want to build this container from?
ARG REPO_NAME=askappy
ARG REPO_VERSION=base

# Record useful metadata: See https://docs.docker.com/config/labels-custom-metadata/
LABEL repo.name="${REPO_NAME}" \
      repo.version="${REPO_VERSION}" \
      repo.maintainer.name="Wasim Raja" \
      repo.maintainer.email="Wasim.Raja@csiro.au" \
      repo.description="The image provides the software env needed for building casacore."

ARG REQFILE=requirements.txt
ARG BUILD_DIR=/_build

ARG MEASURES_DATA="/_build/measures_data"
ARG INSTALL_DIR=/_install

RUN mkdir -p ${INSTALL_DIR} ${BUILD_DIR}
RUN mkdir -p ${MEASURES_DATA}


# MPI:
#ARG MPICH_VERSION="3.4.2"
#ARG CPU_CORE_COUNT=1
#ARG MPICH_CONFIGURE_OPTIONS="--enable-fast=all,O3 --with-device=ch4:ofi --prefix=/usr/local"
#ARG MPICH_MAKE_OPTIONS="-j${CPU_CORE_COUNT}"
#RUN mkdir -p ${BUILD_DIR}/mpich-build \
#    && cd ${BUILD_DIR}/mpich-build \
#    && wget https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
#    && tar xzf mpich-${MPICH_VERSION}.tar.gz \
#    && cd mpich-${MPICH_VERSION}  \
#    && ./configure ${MPICH_CONFIGURE_OPTIONS} \
#    && make ${MPICH_MAKE_OPTIONS} \
#    && make install \
#    && ldconfig \
#    && cd / \
#    && rm -rf ${BUILD_DIR}/mpich-build

## ***************** Build from code base ********************
RUN mkdir /app \
    && cd app \
    && mkdir vast_fast \
    && cd /app/vast_fast

COPY . .

RUN /venv/bin/python3.8 -m pip install --upgrade pip
RUN pip3 install poetry
RUN poetry install
# RUN pip3 install pyregion
## +++++++++++++++++ Build casacore from source +++++++++++++++++++++
## Build casacore first - get measures data as a first step:
## TODO: We need a way to link measures data to the updated version
## kept on host machines.
#RUN cd ${MEASURES_DATA} \
#    && wget ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar \
#    && tar xf *.ztar
#RUN cd ${BUILD_DIR} \
#    && git clone https://github.com/casacore/casacore \
#    && cd casacore \
#    && mkdir build \
#    && cd build \
#    && cmake -DDATA_DIR=${MEASURES_DATA} -DUSE_OPENMP=OFF \
#    -DBUILD_PYTHON3=ON -DUSE_THREADS=ON .. \
#    && make -j4 install
## +++++++++++++++++++ End of casacore build +++++++++++++++++++++++
#
### Fetch the requirement file from the repo and install requirements:
#COPY ${REQFILE} /requirements.txt
#
#RUN pip3 install --no-cache-dir --no-deps -r requirements.txt
#RUN pip3 install --default-timeout=120 --no-binary python-casacore python-casacore
#
## Install modular casatools and caatasks - can be invoked from python environments:
## See: https://casa.nrao.edu/casadocs/casa-6.1.0/usingcasa/obtaining-and-installing
#RUN mkdir -p ${BUILD_DIR}/${REPO_NAME} \
#    && cd ${BUILD_DIR}/${REPO_NAME} \
#    && pip3 install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools \
#    && pip3 install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks \
#    && pip3 install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casampi \
#    && pip3 install casaplotms \
#    && pip3 install casaviewer
#
## Define the measures data area in the root's home area:
#ARG CASA_CONFIG_DIR="/root/.casa"
#RUN mkdir -p ${CASA_CONFIG_DIR}
#COPY config.py ${CASA_CONFIG_DIR}/.
#COPY Dockerfile /Dockerfile
#
## Replicate the policy for ImageMagick from Setonix setup
#ADD ImageMagick-7_on_setonix_policy.xml /etc/ImageMagick-6/policy.xml
#ENV LD_LIBRARY_PATH="/usr/local/lib"
#
#CMD ["/bin/bash"]

