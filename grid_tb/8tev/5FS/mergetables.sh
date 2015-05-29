#!/bin/bash
rm -f NNPDF23_5f.* CT10_5f.* MSTW2008_5f.*
for i in tanbeta_NNPDF23_???GeV tanbeta_NNPDF23_????GeV; do cat $i | grep -v "#" >>NNPDF23_5f.txt; done
for i in tanbeta_CT10_???GeV tanbeta_CT10_????GeV; do cat $i | grep -v "#" >>CT10_5f.txt; done
for i in tanbeta_MSTW2008_???GeV tanbeta_MSTW2008_????GeV; do cat $i | grep -v "#" >>MSTW2008_5f.txt; done
root -q -b readanal.C\(\"NNPDF23_5f\"\)
root -q -b readanal.C\(\"CT10_5f\"\)
root -q -b readanal.C\(\"MSTW2008_5f\"\)
