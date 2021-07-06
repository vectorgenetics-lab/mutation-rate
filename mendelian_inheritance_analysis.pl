############################################ Mendelian inheritance by ILIYAS ####
#			Run command
#	module load vcftools
#	perl Mendelian_Inheritance_Analysis_AgamFbayes.pl
###########################################################################
#!/usr/bin/perl -w
use strict;
use warnings;
use Vcf;
my$filepath="pconcat";
my$vcffile="pconcat_biallelic_shuffle.vcf.gz";			# contains all four potent male parents
print"\n\nFile Path: $filepath\nFile name: $vcffile\n\n";
my$vcf=Vcf->new(file=>"$filepath/$vcffile");
$vcf->parse_header();
my$format_header=$vcf->format_header();
my @samples = $vcf->get_samples;

my$c=0;	#repesenting row number of the vcf;
my($altFrefM1,$refFaltM1,$altFrefM2,$refFaltM2,$altFrefM3,$refFaltM3,,$altFrefM4,$refFaltM4)=(0,0,0,0,0,0,0,0);
my($SitealtFrefM1,$SiterefFaltM1,$SitealtFrefM2,$SiterefFaltM2,$SitealtFrefM3,$SiterefFaltM3,$SitealtFrefM4,$SiterefFaltM4)=(0,0,0,0,0,0,0,0);
my$recdcount=0;
my@ColumnName=();
for(my$i=0;$i<=41;$i++)
{
	my$columnsn=$vcf->get_column_name($i);
	push(@ColumnName,$columnsn);
	
}
print"\nTotal columns in '$vcffile' ".scalar @ColumnName."\n@ColumnName\n\n";
print"\nAll samples ".scalar @samples."\n@samples\n\n";

my$outfprefix=$vcffile;
$outfprefix=~s/\.vcf.*//g;
my$outfile=$outfprefix."_mendelian.vcf";
my%CHRSITES=();
open(OUT,">$outfile") or die "$outfile cant open for writing:$!\n";
print OUT $format_header;
while(my$rowdata=$vcf->next_data_array)		#my$vcfline= $vcf->next_line()
{$c++;
	#print"$c\t$$rowdata[0]\t$$rowdata[1]\t$$rowdata[2]\t$$rowdata[3]\t$$rowdata[4]\t$$rowdata[5]\t$$rowdata[6]\t$$rowdata[8]\t$$rowdata[9]\t$$rowdata[10]\t$$rowdata[11]\n\n\n";
	my@rowrec=@$rowdata;
	my@sampletype=splice(@rowrec,9);
	my($gtfp,$gtmp1,$gtmp2,$gtmp3,$gtmp4)=&parentsgt(\@sampletype,\@samples);
	if($gtfp=~/(1(\/|\|)1)/ and $gtmp1=~/(0(\/|\|)0)/ and $gtmp1!~/($gtmp2|$gtmp3|$gtmp4)/)	#ie genotype description where FP like 1/1 or 1|1 and MP1 like 0/0 or 0|0, but PM2 not like 0/0 or 0|0
	{$SitealtFrefM1++;
		#print "$gtfp  $gtmp1   $gtmp2   $gtmp3   $gtmp4\n";
		#print"$sampletype[0]\t$sampletype[1]\t$sampletype[2]\t$sampletype[3]\t$sampletype[4]\n\n";
		my$nonhetero=offspringtest(\@sampletype,\@samples);	# Tested all offspring are heterozygous for a site (a row)
		if($nonhetero==0)	# one non heterozygous offspring is allowed because #F1foc07f showing 95% nocall
		{$altFrefM1++;$recdcount++;
			my$i=-1;
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			#print "FP: $gtfp\tMP1: $gtmp1\tMP2: $gtmp2\n";
			foreach (@sampletype)	#for printing data of a row to check heterozygous offspring
			{$i++;
				#print "$i\t$samples[$i]\t$_\n";	
			}
			#print"\n\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
		#last if($altFrefM1>=1);
	}
	if($gtfp=~/(0(\/|\|)0)/ and $gtmp1=~/(1(\/|\|)1)/ and $gtmp1!~/($gtmp2|$gtmp3|$gtmp4)/)
	{$SiterefFaltM1++;	
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$refFaltM1++;$recdcount++;
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
	if($gtfp=~/(1(\/|\|)1)/ and $gtmp2=~/(0(\/|\|)0)/ and $gtmp2!~/($gtmp1|$gtmp3|$gtmp4)/)
	{$SitealtFrefM2++;		
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$altFrefM2++;$recdcount++;
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
	if($gtfp=~/(0(\/|\|)0)/ and $gtmp2=~/(1(\/|\|)1)/ and $gtmp2!~/($gtmp1|$gtmp3|$gtmp4)/)
	{$SiterefFaltM2++;
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$refFaltM2++;$recdcount++;	
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
###############################
	if($gtfp=~/(1(\/|\|)1)/ and $gtmp3=~/(0(\/|\|)0)/ and $gtmp3!~/($gtmp1|$gtmp2|$gtmp4)/)
	{$SitealtFrefM3++;
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$altFrefM3++;$recdcount++;	
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
	if($gtfp=~/(0(\/|\|)0)/ and $gtmp3=~/(1(\/|\|)1)/ and $gtmp3!~/($gtmp1|$gtmp2|$gtmp4)/)
	{$SiterefFaltM3++;
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$refFaltM3++;$recdcount++;	
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
	if($gtfp=~/(1(\/|\|)1)/ and $gtmp4=~/(0(\/|\|)0)/ and $gtmp4!~/($gtmp1|$gtmp2|$gtmp3)/)
	{$SitealtFrefM4++;
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$altFrefM4++;$recdcount++;	
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}
	if($gtfp=~/(0(\/|\|)0)/ and $gtmp4=~/(1(\/|\|)1)/ and $gtmp4!~/($gtmp1|$gtmp2|$gtmp3)/)
	{$SiterefFaltM4++;
		my$nonhetero=offspringtest(\@sampletype,\@samples);
		if($nonhetero==0)
		{$refFaltM4++;$recdcount++;	
			my$rec=join("\t",@$rowdata);
			print OUT"$rec\n";
			if(exists $CHRSITES{$$rowdata[0]}){$CHRSITES{$$rowdata[0]}++;}else{$CHRSITES{$$rowdata[0]}=1;}
		}
	}	
}
close OUT;
my$altrefhomocomb=$SitealtFrefM1+$SiterefFaltM1+$SitealtFrefM2+$SiterefFaltM2;
print"All biallelic sites in the input file '$filepath/$vcffile':$c\n\n";
print "Total sites in vcf file as ref and alt homozygous combinations among all parents: $altrefhomocomb\n";
print"\nA file of all Mandilian Inheritence based sites ($recdcount) is written  in the current directory\n\n";

my%STATS=();
my$mendelianchart="Homozygous type\tHomozygous parents\tHeterozygous offspring\n";
$STATS{"Alt FP and Ref MP1 (1/1 X 0/0)\t$SitealtFrefM1"}="$altFrefM1";
$STATS{"Ref FP and Alt MP1 (0/0 X 1/1)\t$SiterefFaltM1"}="$refFaltM1";
$STATS{"Alt FP and Ref MP2 (1/1 X 0/0)\t$SitealtFrefM2"}="$altFrefM2";
$STATS{"Ref FP and Alt MP2 (0/0 X 1/1)\t$SiterefFaltM2"}="$refFaltM2";
$STATS{"Alt FP and Ref MP3 (1/1 X 0/0)\t$SitealtFrefM3"}="$altFrefM3";
$STATS{"Ref FP and Alt MP3 (0/0 X 1/1)\t$SiterefFaltM3"}="$refFaltM3";
$STATS{"Alt FP and Ref MP4 (1/1 X 0/0)\t$SitealtFrefM4"}="$altFrefM4";
$STATS{"Ref FP and Alt MP4 (0/0 X 1/1)\t$SiterefFaltM4"}="$refFaltM4";

foreach my$data(sort {$STATS{$b} <=> $STATS{$a}}keys %STATS)
{
		my($homotype,$homoparents)=split(/\t/,$data);
		my$hetoffspring=$STATS{$data};
		my$coverdpercent=sprintf("%.2f",($hetoffspring*100)/$homoparents);
		$mendelianchart.= "$homotype, $coverdpercent\t$homoparents\t$hetoffspring\n";
		print"No. of sites with all heterozygous offspring are $hetoffspring among $homoparents homozygous '$homotype' parents, that is ".$coverdpercent ." %\n";
		#print"A total $hetoffspring ($coverdpercent %) sites represent all heterozygous offspring among $homoparents sites of '$homotype' homozygous parents\n";
}
my$outstat=$outfprefix.'_MendelianStats.txt';
open(OUT2,">$outstat") or die"file cant written\n";
print"\n\nData for bar diagram\n\n$mendelianchart\n\n";
print OUT2"Data for bar diagram\n\n$mendelianchart\n\n";
print "\nAll Mendilian sites found in different chromosomes\n\nCHR\tNo. of sites\n";
print OUT2"\n\n\nAll Mendilian sites found in different chromosomes\n\nCHR\tNo. of sites\n";
foreach (keys %CHRSITES)
{
	print OUT2"$_\t$CHRSITES{$_}\n";
	print "$_\t$CHRSITES{$_}\n";
}
##################################################################
sub parentsgt()
{
	#GT:AD:DP:GQ:PGT:PID:PL:PS			GATK format	
	#GT:DP:AD:RO:QR:AO:QA:GL			Freebayes format
	my($formatdata,$sample)=@_;
	my@fp=split(/\:/,$$formatdata[0]);
	my@mp1=split(/\:/,$$formatdata[1]);
	my@mp2=split(/\:/,$$formatdata[2]);
	my@mp3=split(/\:/,$$formatdata[3]);
	my@mp4=split(/\:/,$$formatdata[4]);
	return($fp[0],$mp1[0],$mp2[0],$mp3[0],$mp4[0]);
}
sub offspringtest
{
	my($sampdata,$sample)=@_;
	my$GT=0;
	my$i=-1;
	foreach(@$sampdata)
	{$i++;
		if($i>4)
		{
			#GT:AD:DP:GQ:PGT:PID:PL:PS
			my@desc=split('\:',$_);
			my$gt=$desc[0];
			#print"$gt\t$$sampdata[$i]\t$$sample[$i]\n";	
			unless($gt=~/(0(\/|\|)1)/g or $gt=~/(1(\/|\|)0)/g) # ie 0/1 or 0|1 1/0 or 1|0 all non heterozygous check
			{
				$GT++;		
			}
		}
	}
	#print"\n\n";
	return ($GT);
}
