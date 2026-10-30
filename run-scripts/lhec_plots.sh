#!/bin/bash

#./plottool unsub resumbk pdrc z2imp unb 0.075 11.2959 1.225 0.01 41.7653 | tee ./out/lhec-lightq-urpi-x001.dat
#./plottool unsub kcbk pdrc z2imp unb 0.075 48.9819 1.25 0.01 48.6216 | tee ./out/lhec-lightq-ukpi-x001.dat
#./plottool unsub trbk pdrc z2imp unb 0.086 68.9756 1.82 0.01 59.9899 | tee ./out/lhec-lightq-utpi-x001.dat

    #sdrc
./plottool unsub resumbk sdrc z2imp unb 0.079 0.707107 1.8 0.01 29.3386 | tee hi.dat
./plottool unsub kcbk sdrc z2imp unb 0.072 1.97453 1.566 0.01 32.356 | tee ./out/lhec-lightq-uksi-x001.dat
./plottool unsub trbk sdrc z2imp unb 0.028 212.398 1.5 0.01 132.237 | tee ./out/lhec-lightq-utsi-x001.dat

# hera
    # pdrc
./plottool unsub resumbk pdrc z2imp unb 0.0964 1.21086 0.98 0.01 19.6678 | tee ./out/lhec-hera-urpi-x001.dat
./plottool unsub kcbk pdrc z2imp unb 0.085 4.24795 1.0 0.01 26.5772 | tee ./out/lhec-hera-ukpi-x001.dat
./plottool unsub trbk pdrc z2imp unb 0.09 9.75463 1.38 0.01 42.8652 | tee ./out/lhec-hera-utpi-x001.dat

