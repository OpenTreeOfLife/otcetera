sh bootstrap.sh || exit
mkdir buildgcc
(cd buildgcc || exit ; \
 bash ../reconf-gcc-travis.sh && make && make)
