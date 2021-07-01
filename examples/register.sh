#!/bin/bash
mkdir -p maps
../register -pairs pairs.lst -input images/ -tif -output maps/ -distortion 2.0 -output_level 6 -depth 6 -quality 0.5 -summary maps/summary.out -input_map cmaps/
