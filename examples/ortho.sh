#!/bin/bash
mkdir -p aligned_yz
../ortho -input aligned/ -image_list images.lst -output aligned_yz/x -format z%d -memory 6000 -yz_images -x 0-1024
