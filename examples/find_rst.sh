#!/bin/bash
mkdir -p cmaps
../find_rst -pairs pairs.lst -tif -images images/ -output cmaps/ -max_res 1024 -scale 1.0 -tx -50-50 -ty -50-50 -summary cmaps/summary.out
