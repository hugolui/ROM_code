#!/bin/bash
sh clean.sh
python reconst_best_model.py &&
sh compile.sh &&
./reconst.out






