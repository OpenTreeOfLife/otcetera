sh bootstrap.sh || exit
mkdir buildgcc
(cd buildgcc && bash ../reconf-gcc.sh && make -j4 && make)
