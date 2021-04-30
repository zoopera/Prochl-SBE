#!/bin/bash

# this scriptis calculate the PNC(dC) and PNR(dR) values for each gene family in the control and target clade
# the output is ./dRdC_dat/HON*/PNC_PNR.target_*.control_*.[charge|MY].txt

WORK_DIR=$(pwd)
TARGET_CLADES=$(ls clades.target_*.control*.txt | cut -d"." -f2 | cut -d"_" -f2- | sort | uniq)

for TARGET in ${TARGET_CLADES[@]} ; do
  cd ${WORK_DIR}
  CONTROL_CLADES=$(ls clades.target_${TARGET}.control_*.txt | cut -d"." -f3 | cut -d"_" -f2- | sort | uniq)
  for CONTROL in ${CONTROL_CLADES[@]} ; do
    echo ":::::: Target = $TARGET, Control = $CONTROL ::::::"
    while read LINE ; do
      if [[ ${LINE:0:3} == 'HON' ]] ; then

        HON=$(echo $LINE | cut -d"=" -f1)
        HON_DIR=$(echo $LINE | cut -d"=" -f2)
        echo ">>> $HON"
        cd ${WORK_DIR}/dRdC_dat/${HON}

        # calculate average dR and dC for each gene
	COMB="target_${TARGET}.control_${CONTROL}"
        for CLASS in 'charge' 'MY' ; do
          # calculate_average_PNC_PNR.pl <classification charge|MY> <clade file> <outfile>
          ${WORK_DIR}/scripts/calc_avg_PNC_PNR_per_genfam.py  $CLASS  ${WORK_DIR}/clades.${COMB}.txt  PNC_PNR.${COMB}.${CLASS}.txt
        done
	num_gene1=$( wc -l PNC_PNR.${COMB}.charge.txt | cut -d" " -f1 )
	num_gene2=$( wc -l PNC_PNR.${COMB}.MY.txt | cut -d" " -f1 )
	printf "    Number of genes used (charge and MY): %s  %s\n" $num_gene1 $num_gene2
      fi
    done < ${WORK_DIR}/dRdC_pipeline.cfg
  done
done
