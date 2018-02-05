FROM ubuntu:16.04

RUN apt-get update

RUN \
  apt-get install -y g++-5 \
  && apt-get install -y autotools-dev \
  && apt-get install -y autoconf \
  && apt-get install -y cmake \
  && apt-get install -y libtool \
  && apt-get install -y libboost-all-dev \
  && apt-get install -y git-all \
  && apt-get install -y python-setuptools python-dev build-essential \
  && apt-get install -y libssl-dev \
  && apt-get install -y libcurl4-openssl-dev \
  && apt-get install -y libcrypto++-dev


RUN \
  git clone --recursive https://github.com/corvusoft/restbed.git \
  && mkdir restbed/build \
  && cd restbed/build  \
  && cmake -DBUILD_SSL=YES .. \
  && make -j8 \
  && make install \
  && cd -

RUN \
  apt-get -y install python-pip \
  && pip install --upgrade pip \
  && pip install requests

RUN \
  git clone https://github.com/mtholder/otcetera.git \
  && cd otcetera \
  && git checkout -b ws origin/ws \
  && bash bootstrap.sh \
  && mkdir build \
  && cd build \
  && export CPPFLAGS="-I/restbed/distribution/include" \
  && export LDFLAGS="-L/restbed/distribution/library" \
  && bash ../reconf-gcc-docker.sh \
  && make -j8 \
  && make check \
  && make install \
  && cd -

EXPOSE 1984

CMD /usr/local/bin/otc-tol-ws otcetera/data/ex-tax-1 -D otcetera/data/ex-synth-par
