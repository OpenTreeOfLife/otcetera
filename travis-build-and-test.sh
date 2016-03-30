sh bootstrap.sh || exit
mkdir buildgcc
(cd buildgcc || exit ; \
 bash ../reconf-gcc-travis.sh || cat config.log; \
 make -j4 && make)
