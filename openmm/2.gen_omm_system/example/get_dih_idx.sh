#!/bin/bash

###
### read dihedrals_to_scale.names.txt and look up
### the atomic indices from the psf (minus 1 for omm 0-indexed)
### prints dihedrals_to_scale.txt
###
### Usage: ./script psf_filename (defaults to patch.psf)
###

# read in psf file
psffile=$1
if [ "$psffile" == "" ]; then
    echo "Getting atomic index numbers from patch.psf"
    psffile=patch.psf
    if [ ! -e $psffile ]; then echo "$psffile not found"; exit; fi
else
    echo "Getting atomic index numbers from $psffile"
fi
if [ ! -e dihedrals_to_scale.names.txt ]; then
    echo "The dihedrals_to_scale.names.txt file not found"
    exit
fi

# for each line in psffile, grep out the atom index number (minus 1)
dlines=`more dihedrals_to_scale.names.txt | wc -l`
for (( a=1; a <= ${dlines}; a++)); do
    tline=(`head -${a} dihedrals_to_scale.names.txt | tail -1`)
    nline="${tline[0]}"
    for b in {1..4}; do
        atidx=`grep "${tline[$b]}" $psffile | awk '{print($1)}'`
        atidx=`echo "scale=0; $atidx - 1" | bc`
        nline="$nline $atidx"
    done
    #echo $nline
    # save to disk too
    if [ $a == 1 ]; then echo $nline > dihedrals_to_scale.txt
    else echo $nline >> dihedrals_to_scale.txt
    fi
done

