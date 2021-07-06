# Mutation Rate
Reference for the paper:

## Spontaneous mutation rate estimation of the principal malaria vector Anopheles coluzzii and Anopheles stephensi
Iliyas Rashid, Melina Campos, Travis Collier, Marc Crepeau, Allison Weakley, Hans Gripkey, Yoosook Lee, Hanno Schmidt, Gregory C. Lanzaro

We estimated the nuclear mutation rate per generation in the malaria vectors Anopheles coluzzii and Anopheles stephensi by detecting de novo genetic mutations. 
 
### Male parent identification
-	Relatedness (--relatedness2 from VCFtools, http://vcftools.sourceforge.net/):
  
-	Mendelian heritage: 
Perl script `mendelian_inheritance_analysis.pl`

### De novo mutation detection
-	De novo mutation detection 
Perl script `denovo_mutation_detection.pl”` : this is the main program for de novo mutation detection

### Synthetic mutation insertion
- Insert mutation using BAMSurgeon (https://github.com/adamewing/bamsurgeon/)
`run_bamsurgeon.sh`

- Fix tags that are not handled correctly by BAMsurgeon 
`run_fixRG.sh`

- Convert the handled correctly by bamsurgeon. Finally, the ynthetically mutated BAM files back into FASTQ 
`run_bamtofastq.sh`
