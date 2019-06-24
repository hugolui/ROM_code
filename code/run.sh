#!/bin/bash
cp inputs.inp pod/ &&
cd pod/ &&
sh clean.sh &&
sh compileOpt.sh &&
./pod && 
cd .. &&
cd regression/compact_scheme/ &&
sh compile.sh &&
./calc_derivative.out &&
cd .. &&
cd deep_learning/ &&
python3 DNN_regression.py &&
python3 MAE.py &&
cd .. &&
cd .. &&
cd reconst &&
python3 reconst_best_model.py 







