#!/bin/bash

rm -f grid.txt
touch grid.txt

echo XX mhp tb xsec unc- unc+ >>grid.txt
#for i in `seq 200 20 1000`; do
for i in `seq 200 20 1000` 1200 1400; do
  echo ZZZ $i
  root -b plot_tb.C\($i\)<<<"cout << endl;" | grep XX | sed s#"XX "##g  >>grid.txt
done
