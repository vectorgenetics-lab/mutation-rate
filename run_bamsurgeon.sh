#!/bin/bash

module load bamsurgeon
module load samtools

for samp in 'F1bait09m' 'F1foc01f' 'F1foc02f' 'F1foc03f' 'F1foc04f' 'F1foc08m' 'F1foc09m' 'F1foc10m' 'F1foc11m' 'F1foc12m'; do
    echo
    echo "###############################"
    echo $samp
    echo "###############################"
    addsnv.py -p 8 --mindepth 1 --ignoresnps --ignorepileup --force --aligner mem --alignopts M: --picardjar "/software/picard-tools/1.139/static/picard.jar" -r "/share/lanzarolab/archive/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa" -f "../Agam_bams/${samp}.bam" -v "${samp}_minDP10_modeDPdeltaNone.varfile" -o "${samp}_MUT_minDP10_modeDPdeltaNone.bam" 2>&1 | tee "${samp}_MUT_minDP10_modeDPdeltaNone.log"
done


