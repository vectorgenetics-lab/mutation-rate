############# This script is written by ILIYAS for calling candidate mutation ####################################################################################################################################################
#This program can be used after preparation of input file as:
	# Biallelic SNPs 
	# Shuffled Sample in in order parents-> bait -> focal
	# Accessibility.genome mapped
#Run command:		$perl DenovoMutationRefhomo_AgamFbayes_noacc.pl /share/lanzarolab/users/mdelima/mutation-rate/pconcat-order-pass.vcf.gz
##########################################################################################################################################################################################################################
#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw(first);
my$BaitMaxHeteroAllowed=0;	 #change this number to allow number of hetrozygot bait ie. value 0 means no heterozygot bait and 3 means three heterozygot bait is allowed in the analysis
my$vcffile=$ARGV[0];
my@vcfheader=();
my@vcfcolumn=();
my($biallelicsites,$impureparents,$pureparents,$heterobait,$noheterobait,$callinfocal,$nocallinfocal)=(0,0,0,0,0,0,0);
my$infile='';
#-------------------------- Set parameters 	#permissive=27.28   and restrict=31.25
my$FocalAltADpct=27;				# minimum percentage of reads supporting the alternative allele depth(ALTDP)at a focal heterozygous site
my$maximuFocalAltADpct=100-$FocalAltADpct;	# Maximum percentage of reads supporting ALTDP
my$minParentsDP=10;				# minimum DP of either homozygous parent (pure) 
#-------------------------- ----------------#
my$focald=10;
our$count=0;
our$singlefocalmut=0;
our$twofocalmut=0;
our$mutation=0;
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	F0female	F0male1	F1bait01f	F1bait02f	F1bait03f	F1bait04f	F1bait05f	F1bait06f	F1bait07f	F1bait08m	F1bait09m	F1bait10m	F1bait11m	F1bait12m	F1bait13f	F1bait14m	F1bait15m	F1bait16m	F1bait17m	F1foc01f	F1foc02f	F1foc03f	F1foc04f	F1foc05f	F1foc06f	F1foc07f	F1foc08m	F1foc09m	F1foc10m	F1foc11m	F1foc12m	F1foc13m
if(!$vcffile or $vcffile !~/\.vcf\.gz/g)
{
	print STDERR "\n\nPlease enter a zipped vcf file e.g. extension '.vcf.gz'\n\n";
	exit;
}
else{
	print"\n****************************************************** Input file **************************************************************\n";
	my$basename = basename($vcffile);
	my$dirname  = dirname($vcffile);
	$infile=$basename;
	print"Parent directory path: $dirname\n";
	print"\nVCF file name: $basename\n\n";
}
my%Mutations=("3L\t10978191"=>1,"3L\t30182152"=>2,"3L\t32059779"=>3,"3R\t18382262"=>4,"3R\t49470109"=>5,"2L\t21501855"=>6, "2L\t6472372"=>7,"2R\t32051427"=>8);
my$fprefix=$infile;
$fprefix=~s/\.vcf.*//g;
my$outfile=$fprefix.'_MutationAnalysis_'.$BaitMaxHeteroAllowed.'bait.vcf';
open(OUT,">$outfile") or die "File '$outfile' cant open for writing:$!\n"; 
my$refhomoparentsfile=$fprefix.'_RefHomozygousParents.vcf';
open(HOMOPAR,">$refhomoparentsfile") or die "File '$refhomoparentsfile' cant open for writing:$!\n"; 
#open(NOHOMOPAR,">NonRefHomozygousParents.vcf") or die "File 'NonRefHomozygousParents.vcf' cant open for writing:$!\n"; 
my$header='';
my$coltitle='';
my@focalname=();
my$x=0;
open(NEWMUTLIST,">DetectedPointMutationsAg.vcf") or die"File 'DetectedPointMutationsAg.vcf' cant open for writing:$!\n";
open(IN, "gunzip -c $vcffile |") or die "gunzip $vcffile: $!";
while( my$l=<IN> )
{
    chomp $l;
	if($l=~/^\#\#/g)
	{
		$header.="\n$l";
	}
	elsif($l=~/^#C/)
	{
		print OUT"$l\n";
		print HOMOPAR"$l\n";
		#print NOHOMOPAR"$l\n";
		my@cols=split('\t',$l);
		my$size=scalar @cols;
		my($s1,$s2,$s3)=&sampletype(@cols);		
		my$samplesize =scalar @$s1;
		my$baitsamplesize=scalar @$s2;
		@focalname= @$s3;
		my$focalsamplesize=scalar @focalname;	
		print "Number of All samples: $samplesize (Total vcf column: $size)\n@$s1\n";
		print "\nNumber of All bait samples: $baitsamplesize\n@$s2\n\n";
		print "\nNumber of All focal samples: $focalsamplesize\n@$s3\n";		
		print"********************************************************************************************************************************\n\n\n";
		print NEWMUTLIST"$l\n";
		#exit;
	}
	else{$biallelicsites++;
		my@sample=split('\t',$l);
		if(exists $Mutations{"$sample[0]\t$sample[1]"})
		{
			print NEWMUTLIST"$l\n";
		}
		
		my$parentcheck=parentcheck($sample[9],$sample[10],$minParentsDP);
		if($parentcheck > 0)	#true parents check #filtered both parents on applied threshold (i) pure homozygous(ii) min depth >10 and (iii) gene quality > 20
		{$pureparents++;
			print HOMOPAR"$l\n";
			my$length=scalar @sample;
			#print"$length sample size\n";
			#exit;
			my($all,$bait,$focal)=&sampletype(@sample);
			my@allsample=@$all;
			my@baitsample=@$bait;
			my@focalsample=@$focal;
			my$heterozygot=baitcheck(\@baitsample,\$BaitMaxHeteroAllowed);	# unpahsed heterozygosity (GT:0/1) testing in the bait samples
			my$flag=0;
			my$atleast=0;
			my$focalheterosamp='';
			if($heterozygot==0)	# if heterozygosity is false in all bait sample 
			{$noheterobait++;
				if($sample[0] eq 'X')
				{
					#print "$x\t$l\n\n\n";
					my$focalcall=ChromosomeX(\@focalsample,\@focalname,\@sample);
					if($focalcall>0){$callinfocal++; print OUT"$l\n";}
				}
				elsif($sample[0] eq 'UNKN')
				{
					next;
				}
				else{
					#next;
					my$nocal=focalcheck(\@focalsample,\'./.',\'1/1' ,\'1|1');
					if($nocal==0) # no call is false
					{$callinfocal++;
						my$indx=0;
						foreach my$s(@focalsample)
						{	my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($s);
							if($gt eq '0/1' or $gt eq '0|1') 	# heterozygosity is true with altleast 50% AD in atleast one bait sample for an allele in vcf
							{$atleast++;
								#my$indx = first { $focalsample[$_] eq $s } 0..$#focalsample;
								my$sampname=$focalname[$indx];
								$focalheterosamp.="$s<-$sampname\t";
								$flag=1;							
							}
							$indx++;
						}			
						if($flag==1 && $atleast <=2)
						{$mutation++;
							
							#print OUT"$l\n";
							$focalheterosamp=~s/(^.+)\t$/$1/g;
							my@focsamp=split(/\t/,$focalheterosamp);
							my$heterosize=scalar @focsamp;
							if($heterosize==1)
							{$singlefocalmut++;
								my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($focsamp[0]);
								if($dp>=$focald && $adalt>= $dp*($FocalAltADpct/100) && $adalt <= $dp*($maximuFocalAltADpct/100))
								{$count++;
									print OUT"$l\n";
									print "$count\t$sample[0]\t$sample[1]\t$sample[2]\t$sample[3]\t$sample[4]\t$focalheterosamp\t$heterosize\n";						
								}							
							} 
							else{$twofocalmut++;						
								my($gt1,$adref1,$adalt1,$dp1,$gq1,$pl1)=sampledesc($focsamp[0]);
								my($gt2,$adref2,$adalt2,$dp2,$gq2,$pl2)=sampledesc($focsamp[1]);
								if(($dp1>=$focald && $adalt1>= $dp1*($FocalAltADpct/100) && $adalt1 <= $dp1*($maximuFocalAltADpct/100)) or ( $dp2>=$focald && $adalt2 >=$dp2*($FocalAltADpct/100) && $adalt2 <= $dp2*($maximuFocalAltADpct/100)))
								{$count++;
									print OUT"$l\n";
									print "$count\t$sample[0]\t$sample[1]\t$sample[2]\t$sample[3]\t$sample[4]\t$focalheterosamp\t$heterosize\n";						
								}						
							}
							$flag=0;
							$atleast=0;
						}
					}	#close focal test
				else{$nocallinfocal++;}
				} #close esle of non x chromosome
			}	#close bait check
			else{$heterobait++;}
		}	#parent check closed
		else{$impureparents++;
			#print NOHOMOPAR"$l\n";
		}	
	}
}
#print"$header\n$coltitle\n";
if($BaitMaxHeteroAllowed==0){$BaitMaxHeteroAllowed ="Not allowed any heterozygous bait";}
elsif($BaitMaxHeteroAllowed==1){$BaitMaxHeteroAllowed="$BaitMaxHeteroAllowed heterozygous bait individual allowed"}
else{$BaitMaxHeteroAllowed="$BaitMaxHeteroAllowed heterozygous bait individuals are allowed"}
print "\n\n----------------------------------------------------- Processing summary -------------------------------------------------------\n";
print "$biallelicsites are biallelic sites in the input vcf file '$infile'\n";
print "$pureparents are ref-homozygous parents (0/0 or 0|0) on minimum DP>=$minParentsDP, ready to bait test\n";
print "$noheterobait are filtered sites without non-reference reads in any bait for the focal test.\n";
print "$callinfocal sites are showing GATK call (0/0, 0|0, 0/1 and 0|1) in focal sample after focal test\n";
print "$mutation are candidate mutation sites among $callinfocal GATK call but maximum two heterozygote focals are allowed.\n";
print "$singlefocalmut sits showing a mutation in a single focal sample and $twofocalmut sits showing a mutation in two focal samples of the total $mutation.\n";
print "$count filtered mutation sites using 'AAF>=$FocalAltADpct%' of the read depth.\n";
print "$BaitMaxHeteroAllowed. Here each allowed heterozygous bait possess a maximum three reads of non-reference allele (ALTAD<4).\n";
print "-------------------------------------------------------------------------------------------------------------------------------\n";
print"\nFiltered based on:\n\tMinimum read depth of homozygous parents: $minParentsDP\n\tPecentage range of Alternate Allele Frequency (AFF)>=$FocalAltADpct% of DPHF && AFF=<$maximuFocalAltADpct% of DPHF\n\n";
print"\nAn output file '$outfile' was written in current working directory\n\n";
print "================================================================================================================================\n\n";
####################### methods #############################################
sub sampletype{
		my@cols=@_;
		my$size=scalar @cols;
		my@baitsamplename=@cols;
		my@focalsamplename=@cols;
		my@samplename=splice(@cols,9,$size);
		my@baitsample=splice(@baitsamplename,11,19);
		my@focalsample=splice(@focalsamplename,30,10);			
		return (\@samplename,\@baitsample,\@focalsample);
}
sub sampledesc
{
	my$sample=shift;
	my($GT,$ADref,$ADalt,$DP,$RO,$QR)=('','','','','');
	#GT:AD:DP:GQ:PL	0/0:50,0:50:77:0,77,1903
	#GT:DP:AD:RO:QR:AO:QA:GL 0/1:15:13,2:13:533:2:78:-2.87328,0,-43.8004
	if($sample=~/^((\d|\.)(\/|\|)(\d|\.))\:(\d+)\:(\d+)\,(\d+)\:(\.|\d+)\:(.+)/g)
	{
		$GT=$1;
		$DP=$5;
		$ADref=$6;
		$ADalt=$7;		
		$RO=$8;
		$QR=$9;		
	}
	#my($GT,$DP,$AD,$RO,$QR,$AO,$QA,$GL)=split(/\:/,$sample);
	#$GT='0/0' if ($GT eq '0|1');
	#my($ADref,$ADalt)=split(/\,/,$AD);
	#print"$sample\t...\tGT:$GT\tAD:$ADref...$ADalt\tDP:$DP\tRO:$RO\tQR:$QR\n";
	return($GT,$ADref,$ADalt,$DP,$RO,$QR);
}
sub parentcheck{
	my($femalesamp,$malesamp,$minPDP)=@_;
	my($fgt,$fadref,$fadalt,$fdp,$fgq,$fpl)=sampledesc($femalesamp); #Female sample
	my($mgt,$madref,$madalt,$mdp,$mgq,$mpl)=sampledesc($malesamp); #Male sample	my$purehomo=0;
	my$purehomo=0;
	if(($fgt eq '0/0' or $fgt eq '0|0')  and ($mgt eq '0/0' or $mgt eq '0|0'))	 #&& $fadalt==0 && $madalt==0) #Allowed only pure homozygous parents ie. GT= 0/0 and AD =n,0:DP=n where n is total number of reads
	{
		if($fgq=~/\d+/ && $mgq=~/\d+/ && $fadalt==0 && $madalt==0)	# test gene quality is a number
		{	if($fdp>=$minPDP && $mdp>=$minPDP)	#use cut off of minDep>=10 and gene quality >=20 for both parents
			{
				#print"\n$femalesamp\n$malesamp\nF: GT:$fgt ADALT:$fadalt  GQ:$fgq\tM: GT:$mgt ADALT:$madalt GQ:$mgq\n";
				$purehomo++;								
			}
		}
		else{	#print "Gene Quality not a number haence excluded:$femalesamp\t$malesamp\n";
		}
	}		
	return ($purehomo);
}
sub baitcheck{
		my($arg1,$arg2)=@_;
		my@bait=@$arg1;
		my$MaxHeteroAllowed=$$arg2;
		#print "@bait\tMaxHeteroAllowed:$MaxHeteroAllowed\n";
		my$size=scalar @bait;
		my$hetero=0;
		my$flag=0;
		my$i=0;
		my($greatadalt,$lessadalt)=(0,0);				
		foreach my$sample(@bait)
		{$i++;
			my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($sample);
			if($gt eq '0/1' or $gt eq '0|1')	# check no. of heterozygous bait
			{
				if($adalt >3)					# count no. of heterozygous bait in which more than 3 nonreference allelic reads (Alt-AD) 
				{
					$greatadalt++;
				}
				else{							# count no. of heterozygous bait in which leass than 3 nonreference allelic reads (Alt-AD)
					$lessadalt++;
				}
			}
			if($gt eq '0/0' or $gt eq '0|0')
			{
				if($adalt>0)
				{
					$hetero++;
				}		
			}
			if($gt eq '1/1' or $gt eq '1|1')
			{
				$hetero++;
			}
		}		
		if($greatadalt>0 || $lessadalt>$MaxHeteroAllowed) # discarded 
		{
			$hetero++;
		}
		if($greatadalt==0 && $lessadalt <=$MaxHeteroAllowed && $hetero==0) #all homozygous bait except maximum allowed bait individual contains less than three nonreference reads 
		{
			$hetero=0;
		}
		#print"\nsample size $size\tTotal hetero zygous unphased genotype:$hetero\n";
		return ($hetero);
}
sub focalcheck{
		my($samp,$no,$unphs,$phs)=@_; 	#\'./.',\'1/1' ,\'1|1'
		my@focal=@$samp;
		my($nocallgt, $unphasegt,$phasegt)=($$no,$$unphs,$$phs);
		my$size=scalar @focal;
		my$nocall=0;
		foreach my$sample(@focal)
		{
			if($sample=~/\.\:\./) #.:.:.:.:.:.:.:.
			{
				$nocall++;
			}
			else{
				my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($sample);
				if(($gt eq $nocallgt) or ($gt eq $unphasegt) or ($gt eq $phasegt))
				{
					$nocall++;
				}
				if($gt eq ('0/0') or $gt eq ('0|0'))
				{
					if($adalt>0)
					{
						$nocall++;
					}
				}
			}
		}
		#print"\nsample size $size\tTotal hetero zygous unphased genotype:$hetero\n";
		return ($nocall);
}
sub ChromosomeX
{
	my($sm1,$sm2,$sm3)=@_;
	my@fclsample=@$sm1;
	my@fclname=@$sm2;
	my@sample=@$sm3;
	my%allfocal=();
	my$i=0;
	my$atleast=0;
	my$xmutation=0;
	my$nocall=0;
	my$flag=0;
	my$focalheterosamp='';	
	foreach(@fclsample)
	{
		my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($_);
		if($fclname[$i]=~/^.+m$/gi)
		{
			if($gt eq './.' or $gt eq '0/1' or $gt eq '0|1')
			{
				$nocall++;
			}			
		}
		else{
			if($gt eq './.' or $gt eq '1/1' or $gt eq '1|1')
			{
				$nocall++;
			}
		}
		$i++;
	}
	#print"\n";
	my$j=0;	
	my$singlemut=0;
	my$twomut=0;
	
	if($nocall==0)
	{
		foreach my$sm(@fclsample)
		{
			my($gt,$adref,$adalt,$dp,$gq,$pl)=sampledesc($sm);
			#print"$j\t$fclname[$j]\t$gt,$adref,$adalt,$dp,$gq,$pl\n";
			my$sampname=$fclname[$j];
			if($sampname=~/^.+\df$/) #for female sample
			{
				if($gt eq '0/1' or $gt eq '0|1')
				{
					if($adalt>= $dp*(30/100) && $adalt <= $dp*($maximuFocalAltADpct/100))	#$adalt>= $dp*($FocalAltADpct/100) && $adalt <= $dp*($maximuFocalAltADpct/100
					{
						$atleast++;
						$focalheterosamp.="$sm<-$sampname\t";
						$flag=1;
					}
				}
							
			}
			else
			{
				if($gt eq '1/1' or $gt eq '1|1' )
				{
					if($dp>3)
					{
						$atleast++;
						$focalheterosamp.="$sm<-$sampname\t";
						$flag=1;
					}
				}
			}
			$j++;
		}
		#print"\n";
		if($flag==1 && $atleast <=2)
		{$xmutation++;$mutation++;
			
			#print OUT"$l\n";
			$focalheterosamp=~s/(^.+)\t$/$1/g;
			my@focsamp=split(/\t/,$focalheterosamp);
			my$heterosize=scalar @focsamp;
			if($heterosize==1)
			{$singlefocalmut++;
				$count++;
				print "$count\t$sample[0]\t$sample[1]\t$sample[2]\t$sample[3]\t$sample[4]\t$focalheterosamp\t$heterosize\n";
			} 
			else{$twofocalmut++;						
				$count++;
				print "$count\t$sample[0]\t$sample[1]\t$sample[2]\t$sample[3]\t$sample[4]\t$focalheterosamp\t$heterosize\n";
			}
			$flag=0;
			$atleast=0;
		}
	}
	#print"$twomut sites showing mutation in two focal samples and $singlemut sites showing mutation in single focal sample of the total $xmutation\n";
	return($xmutation);
}
system("perl meandepth_calculation_mutatedsite_AgamFbayes.pl $outfile");
