sh bootstrap.sh || exit
mkdir buildgcc
(cd buildgcc && bash ../reconf-gcc-travis.sh && make -j4 && make)
