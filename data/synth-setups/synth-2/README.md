I'm not completely sure where the taxonomy came from here.

The synth tree was build using:

$ cd propinquity
$ export OTT_DIR=$HOME/Devel/OpenTree/otcetera/source/data/synth-setups/taxonomy/
$ export PROPINQUITY_OUT_DIR=$HOME/Devel/OpenTree/otcetera/source/data/synth-setups/synth-par/synth-1
$ make clean && ./build-from-newicks.sh ~/Devel/OpenTree/otcetera/source/data/synth-setups/phyloinputs/ranking.txt
$ make all
$ make extra

