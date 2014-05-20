#!/bin/bash
echo $1
echo $2

openfile=$1
numfiles=$2
for (( m=1; m<$numfiles; m++)) ; do
   ../../../utilities/PostProcessing/FldToVtk-g $openfile\.xml $openfile\_$m\.chk
done
