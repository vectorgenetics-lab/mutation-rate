#!/bin/bash

module load samtools
module load picard-tools

for samp in 'F1bait09m' 'F1foc01f' 'F1foc02f' 'F1foc03f' 'F1foc04f' 'F1foc08m' 'F1foc09m' 'F1foc10m' 'F1foc11m' 'F1foc12m'; do
    echo
    echo "###############################"
    echo $samp
    echo "###############################"
samtools sort -m 4G -O bam ${samp}_MUT_minDP10_modeDPdeltaNone.bam |
picard AddOrReplaceReadGroups \
-INPUT /dev/stdin \
-OUTPUT ${samp}_MUT_minDP10_modeDPdeltaNone_fixRG.bam \
-RGID ${samp} \
-RGLB ${samp} \
-RGPL illumina \
-RGPU UNKN \
-RGSM ${samp} \
-SORT_ORDER coordinate \
-VALIDATION_STRINGENCY LENIENT \
-TMP_DIR /tmp/
samtools index -b ${samp}_MUT_minDP10_modeDPdeltaNone_fixRG.bam

done

