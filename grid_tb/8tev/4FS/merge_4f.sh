#!/bin/bash

function transf()
{
  rm -f ${2}
#  for i in ${1}*tanbscan*; do 
  for i in ${1}*tanbscan_mh_???_* ${1}*tanbscan_mh_????_*  ; do 
#    x=
    mass=`echo $i | sed s#${1}_nf4_tanbscan_mh_##g | sed s#_8tev.data##g | sed s#_14tev.data##g | sed s#_13tev.data##g `
    exec<$i
    while read line; do
      tb=`echo $line | awk '{print $1}'`
      xs=`echo $line | awk '{print $2}'`
      xs_pb=$(awk 'BEGIN {print '"$xs/1000.0"'}' < /dev/null)
#      xs_pb=$(echo "scale=6; $xs/1000." | bc -q 2>/dev/null)
      echo $mass $tb  $xs_pb >>${2}
    done
  #  cat $i | sed s#^#$mass'   #'g  >>CT10.txt; 
  done
}

transf ct10 CT10_4f.txt
transf mstw MSTW2008_4f.txt
transf nnpdf NNPDF23_4f.txt

root -l readanal.C\(\"NNPDF23_4f\"\)
root -l readanal.C\(\"CT10_4f\"\)
root -l readanal.C\(\"MSTW2008_4f\"\)
