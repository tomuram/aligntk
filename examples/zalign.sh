#!/bin/bash
mkdir -p amaps
mkdir -p grids
../align -images images/ -image_list images.lst -map_list pairs.lst -maps maps/ -output amaps/ -schedule schedule.lst -incremental -output_grid grids/ -grid_size 1800x1400 -fixed z08
