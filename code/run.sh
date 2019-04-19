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
python DNN_regression.py &&
python MAE.py &&
cd .. &&
cd .. &&
cd reconst &&
python reconst_best_model.py 







