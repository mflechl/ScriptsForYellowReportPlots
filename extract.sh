#!/bin/bash

ecm=13TeV
echo "${ecm}"

echo "#############################################################"
echo "#############################################################"
echo "combined"
echo "#############################################################"
echo "NLO"
echo "#############################################################"
grep "  combined  " 5F_combined_results_${ecm}.txt | grep -v WARNING                                                              | awk {'print $5","'}
echo "#############################################################"
echo "combined total uncertainty -"
echo "#############################################################"
grep -v "  combined  " 5F_combined_results_${ecm}.txt | grep -v "WARNING\|X 0\|available\|Using\|ATTENTION\|error\|Gaussian"      | awk {'print $2","'}
echo "#############################################################"
echo "combined total uncertainty +"
echo "#############################################################"
grep "  combined  " 5F_combined_results_${ecm}.txt | grep -v WARNING                                                              | awk {'print $6","'}
echo "#############################################################"




echo "#############################################################"
echo "#############################################################"
echo "CT10"
echo "#############################################################"
echo "CT10 NLO"
echo "#############################################################"
grep " CT10 " 5F_CT10_results_${ecm}.txt | awk {'print $4","'}
echo "#############################################################"
echo "CT10 total uncertainty -"
echo "#############################################################"
grep -v CT10 5F_CT10_results_${ecm}.txt | awk {'print $4","'}
echo "#############################################################"
echo "CT10 total uncertainty +"
echo "#############################################################"
grep " CT10 " 5F_CT10_results_${ecm}.txt | awk {'print $8","'}
echo "#############################################################"




echo "#############################################################"
echo "#############################################################"
echo "MSTW"
echo "#############################################################"
echo "MSTW NLO"
echo "#############################################################"
grep " MSTW" 5F_MSTW_results_${ecm}.txt | awk {'print $5","'}
echo "#############################################################"
echo "MSTW total uncertainty -"
echo "#############################################################"
grep -v MSTW 5F_MSTW_results_${ecm}.txt | grep -v "WARNING\|X 0\|available\|Using\|ATTENTION\|error\|Gaussian" | awk {'print $6","'}
echo "#############################################################"
echo "MSTW total uncertainty +"
echo "#############################################################"
grep " MSTW" 5F_MSTW_results_${ecm}.txt | awk {'print $10","'}
echo "#############################################################"

echo "#############################################################"
echo "#############################################################"
echo "NNPDF"
echo "#############################################################"
echo "NNPDF NLO"
echo "#############################################################"
grep " NNP" 5F_NNPDF_results_${ecm}.txt | grep -v "WARNING\|X 0\|available\|Using\|ATTENTION\|error\|Gaussian" | awk {'print $5","'}
echo "#############################################################"
echo "NNPDF total uncertainty -"
echo "#############################################################"
grep -v NNP 5F_NNPDF_results_${ecm}.txt | grep -v "WARNING\|X 0\|available\|Using\|ATTENTION\|error\|Gaussian" | awk {'print $6","'}
echo "#############################################################"
echo "NNPDF total uncertainty +"
echo "#############################################################"
grep " NNP" 5F_NNPDF_results_${ecm}.txt | grep -v "WARNING\|X 0\|available\|Using\|ATTENTION\|error\|Gaussian" | awk {'print $10","'}
echo "#############################################################"


echo "#############################################################"
echo "#############################################################"
echo "scale"
echo "#############################################################"
cat 5F_scale_var_new_results_${ecm}.txt | awk {'print $4", "$5","'}
echo "#############################################################"


#cat ct10_8TeV_4F.data | awk '{print $1" "$2", "$3", "$4","}'
#cat mstw08_8TeV_4F.data | awk '{print $1" "$2", "$10", "$11","}' mass cent err+ err-
#cat nnpdf23_8TeV_4F.data | awk '{print $1" "$2", "$9", "$10","}' 
