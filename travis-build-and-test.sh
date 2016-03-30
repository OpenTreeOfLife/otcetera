sh bootstrap.sh || exit
mkdir buildclang
(cd buildclang && bash ../reconf-clang.sh && make -j4 && make)
