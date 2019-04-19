#!/bin/bash

rm nohup.out
nohup ./pod &
tail -f nohup.out

