autoreconf || exit
git clone https://github.com/miloyip/rapidjson.git
mkdir buildclang
(export RAPID_JSON_INC="$PWD/rapidjson/include" ; cd buildclang && bash ../reconf-clang.sh && make -j4 && make)
