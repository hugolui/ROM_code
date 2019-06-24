#!/bin/bash
sh clean.sh
python3 reconst_best_model.py &&
sh compile.sh &&
./reconst.out






