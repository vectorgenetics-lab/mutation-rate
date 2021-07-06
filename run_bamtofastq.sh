#!/bin/bash
# This script will convert sorted BAM files in the current working directory to paired-end FASTQ files.

module load samtools

OUTDIR=/share/lanzarolab/seq/trim/Agam_mute_rate_redux_synthetic

mkdir -p $OUTDIR
for INPUT in *_fixRG.bam ; do

	BASENAME=${INPUT%_MUT_minDP10_modeDPdeltaNone_fixRG.bam}
	mkdir -p $OUTDIR/${BASENAME}
	mkdir -p $OUTDIR/${BASENAME}/run_synthetic
	echo "sorting ${INPUT}"
	samtools sort -@ 8 -n ${INPUT} -o ${INPUT%.bam}_sorted.bam
	echo "converting ${INPUT%.bam}_sorted.bam to fastq"
	samtools fastq -@ 8 ${INPUT%.bam}_sorted.bam \
		-1 $OUTDIR/${BASENAME}/run_synthetic/trim_pe1.fastq.gz \
		-2 $OUTDIR/${BASENAME}/run_synthetic/trim_pe2.fastq.gz \
		-0 $OUTDIR/${BASENAME}/run_synthetic/trim_sing1.fastq.gz \
		-s $OUTDIR/${BASENAME}/run_synthetic/trim_sing2.fastq.gz -n
	echo "@RG\tID:${BASENAME}_run-synthetic\tLB:KHP-1.1-PE150-${BASENAME}\tPL:ILLUMINA\tSM:${BASENAME}\tPU:unkn" > $OUTDIR/${BASENAME}/run_synthetic/_RG_TAGS
	touch $OUTDIR/${BASENAME}/run_synthetic/trim.done

done
