#!/bin/bash
mkdir -p aligned
../apply_map -image_list images.lst -images images/ -maps amaps/ -output aligned/ -memory 1000
