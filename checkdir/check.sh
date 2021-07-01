#!/bin/bash
echo "Checking AlignTK executables for correctness..."

errors=0
rm -rf cmaps maps amaps grids aligned aligned_yz compare.out


echo ""
echo "Checking find_rst"
echo "1000000.0" >compare.out
mkdir -p cmaps
../find_rst -pairs pairs.lst -tif -images ../examples/images/ -output cmaps/ -max_res 1024 -scale 1.0 -tx -50-50 -ty -50-50 -summary cmaps/summary.out
../compare_maps -map1 cmaps/z00_z01.map -map2 cmaps.correct/z00_z01.map -output compare.out
export result=`cat compare.out`
if (( $(echo "$result < 1.0" | bc -l ) )); then
    echo "find_rst passed check."
else
    echo "find_rst failed check."
    errors=1
fi


echo ""
echo "Checking register"
echo "1000000.0" >compare.out
mkdir -p maps
../register -pairs pairs.lst -images ../examples/images/ -tif -output maps/ -distortion 2.0 -output_level 6 -depth 6 -quality 0.5 -summary maps/summary.out -initial_map cmaps.correct/
../compare_maps -map1 maps/z00_z01.map -map2 maps.correct/z00_z01.map -output compare.out
export result=`cat compare.out`
if (( $(echo "$result < 1.0" | bc -l ) )); then
    echo "register passed check."
else
    echo "register failed check."
    errors=1
fi

echo ""
echo "Checking align"
echo "1000000" >compare.out
mkdir -p amaps
mkdir -p grids
../align -images ../examples/images/ -image_list images.lst -map_list maps.lst -maps maps.correct/ -output amaps/ -schedule schedule.lst -incremental -output_grid grids/ -grid_size 1800x1400 -fixed z08 -font ../font.pgm
../compare_maps -map1 amaps/z00.map -map2 amaps.correct/z00.map -output compare.out
result=`cat compare.out`
if (( $(echo "$result < 1.0" | bc -l ) )); then
    echo "align passed check."
else
    echo "align failed check."
    errors=1
fi

echo ""
echo "Checking apply_map"
echo "1000000" >compare.out
mkdir -p aligned
../apply_map -image_list images.lst -images ../examples/images/ -maps amaps.correct/ -output aligned/ -memory 1000 -font ../font.pgm
../compare_images -image1 aligned/z00.tif -image2 aligned.correct/z00.tif -output compare.out
result=`cat compare.out`
if (( $(echo "$result < 1.0" | bc -l ) )); then
    echo "apply_map passed check."
else
    echo "apply_map failed check."
    errors=1
fi

echo ""
echo "Checking ortho"
echo "1000000" >compare.out
mkdir -p aligned_yz
../ortho -input aligned.correct/ -image_list images.lst -output aligned_yz/x -format z%d -memory 1000 -yz_images -x 512
../compare_images -image1 aligned_yz/x512.tif -image2 aligned_yz.correct/x512.tif -output compare.out
result=`cat compare.out`
if (( $(echo "$result < 1.0" | bc -l ) )); then
    echo "ortho passed check."
else
    echo "ortho failed check."
    errors=1
fi

echo ""
if (($errors > 0)); then
   echo "Errors detected in AlignTK"
   exit 1
else
    echo "All checks completed successfully."
fi
exit 0
