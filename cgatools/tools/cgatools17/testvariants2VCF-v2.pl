#!/usr/bin/perl -w
use strict;

#Converts a cgatools testvariant output to a multi-sample VCF file (VCF spec 4.1)
#Requires cgatools to access reference genome encoded in .crr format (hg18.crr or hg19.crr)
#make sure cgatools is in the $PATH
#
#Example testvariant file:
#variantId       chromosome      begin   end     varType reference       alleleSeq       xRef    GS19238 GS19239	GS19240
#6874944 chr5    20584031        20584032        snp     C       T       dbsnp.119:rs10037487	00	01		11
#6874945 chr5    20584031        20584032        sub     C       TA              		01      01      	00
#6874946 chr5    20584031        20584032        sub     C       TG              		00      01      	00
#After converting to VCF:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	GS19238 GS19239        GS19240
#chr5    20584031        6874944;6874945;6874946 AC      AT,ATA,ATG      .       PASS    NS=51;RSID=dbsnp.119:rs10037487,.,.;AF=0.57,0.02,0.00;DB    GT      0/2	./.	1/1
#
#Variants that share the same location (chr,begin,end) will be merged into one locus and their flags (0,1,N) will be converted into genotype calls
#Samples that are positive for more than two alleles within the same locus will be flagged and their genotype calls set to unknown (./.)
#Look at the sample GS19239 as such an example
#
#For non-SNP locus, VCF requires an extra reference base immediately upstream of the variant locus be included in the REF and ALT columns 
#


die "Usage: $0 testvarOutput.txt hg19.crr > vcf.txt 2> runlog.txt\n\t" . 
	"testvarOutput.txt = output file from cgatools testvariants\n\t" . 
	"hg19.crr = reference genome in .crr format (must be the same as used in testvariants)\n\t" . 
	"vcf.txt = converted file in vcf format\n\t" . 
	"runlog.txt = log file\n" unless ($#ARGV == 1);
die "Fail to find input file " . $ARGV[0] . "\n" unless (-e $ARGV[0]);
die "Fail to find .crr file " . $ARGV[1] . "\n" unless (-e $ARGV[1]);
die "Fail to execute cgatools: $!\n" if (system ('cgatools16 > /dev/null 2>&1') != 0);

my (undef, undef, undef, $mday, $mon, $year) = localtime; 
my($timestamp)=sprintf ('%s%02d%02d', $year+1900, $mon+1, $mday);
print <<EOF;
##fileformat=VCFv4.1
##fileDate=$timestamp
##source=CGA Tools v1.6 listvariants/testvariants
##reference=$ARGV[1]
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples with Fully Called Data">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=RSID,Number=A,Type=String,Description="dbSNP rs ids">
##FILTER=<ID=s50,Description="Less than 50% of samples are fully called">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF

&testvar2VCF (@ARGV);

sub testvar2VCF {
#variantId=0       chromosome=1      begin=2   end=3     varType=4 reference=5       alleleSeq=6       xRef=7    GS12885-1100-37-ASM=8
	my ($inFile, $crrFile) = @_;
	my (@samples, %seen);
	open (IN, "$inFile") or die "Fail to open input file $inFile\n";
	while (<IN>) {
		chomp;
		my (@fs) = split (/\t/);
		if (/^variantId/) {
			@samples = @fs[8..$#fs];
			print join ("\t", ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', @samples)), "\n";
			next;
		}
		next unless (/^\d/);
		my ($key) = $fs[1] . '-' . $fs[2] . '-' . $fs[3];
		if (!%seen || exists $seen{$key}) {
			push (@{$seen{$key}}, [@fs]);
			next;
		}
		#has a new locus, need to process saved locus
		&doSavedLocus (\%seen, \@samples, $crrFile);
		#reinitialize with the new locus
		%seen = ();
		push (@{$seen{$key}}, [@fs]);
	}
	close IN;
	&doSavedLocus (\%seen, \@samples, $crrFile); #do the last locus
}

sub doSavedLocus {
	my ($seen, $samples, $crrFile) = @_;
	my ($lines) = values %$seen; #get all the lines that belong to the same locus
	my ($chr, $begin, $end)  = ($lines->[0]->[1], $lines->[0]->[2], $lines->[0]->[3]);
	my ($info_NS) = scalar @$samples; #num. of samples fully called
	die "no_samples,$chr-$begin-$end\n" unless ($info_NS > 0);
	my ($info_DB) = ''; #dbSNP memebership
	my (@info_RSid) = (); #dbSNP rs ids
	my (%ref, @alt, @id, %alleleCalls);
	my ($format) = 'GT'; #GT (genotype) is the only value currently populated for each sample
	my ($qual, $filter, $i) = ('.', 'PASS', 0); #no quality scores provided in the testvar output
	my (@info_AF);
	#initialize every allele count to 0
	for ($i=0; $i<=scalar @$lines; $i++) {$info_AF[$i] = 0;}
	my ($pos) = $begin;
	#check the vartypes
	my (%varTypes) = map {$_->[4] => 1} @$lines; #column 4 in the testvar output is varType
	my (@varTypeKeys) = keys %varTypes;
	#check if this locus is SNP only
	my ($snpOnly) = 0;
	$snpOnly = 1 if ($#varTypeKeys == 0 && $varTypeKeys[0] eq 'snp');
	foreach my $line (@$lines) {
		my ($refBase) = $line->[5];
		my ($altBase) = $line->[6];
		if (!$snpOnly) { #not a SNP only locus
			my ($decodecrrCmd) = join ('', ("cgatools decodecrr --reference $crrFile --range $chr:", $begin - 1 , "-$begin 2> /dev/null"));
			my ($re) = `$decodecrrCmd`;
			die "cgatools decodecrr failed $decodecrrCmd ", join ("\t", @$line), "\n" unless (defined $re);
			chomp ($re);
			$refBase = $re . $refBase;
			$altBase = $re . $altBase;
			$pos = $begin;
		}
		$ref{$refBase} = 1;
		push (@alt, $altBase);
		push (@id, $line->[0]);
		if ($line->[7] =~ /dbsnp/) {
			$info_DB = 'DB';
			$line->[7] =~ s/;/-/g; #testvar use ; as delimiter which is reserved by the INFO column
			push (@info_RSid, $line->[7]);
		} else { 
			push (@info_RSid, '.');
		}
		#collect testvar's allele flags (00, 01, 11, 0N, NN, ...), which starts at column 8.
		for ($i=0; $i<scalar @$samples; $i++) {
			push (@{$alleleCalls{$samples->[$i]}}, $line->[$i+8]);
		}
	}
	$pos = $begin + 1 if ($snpOnly);
	my (@refBases) = keys %ref;
	die "mismatched_ref_base, " . join (",", (@id, @refBases)) . "\n" if ($#refBases != 0);
	my ($idstr) = join (";", @id);
	my ($altstr) = join (",", @alt);
	print join ("\t", ($chr, $pos, $idstr, $refBases[0], $altstr,$qual)), "\t";
	my (@allGTcalls) = ();
	#convert allele flags to genotype calls
	foreach my $sample  (@$samples) {
		my ($sampleAlleleCalls) = $alleleCalls{$sample};
		my ($sampleGTcall);
		if ($sampleAlleleCalls->[0] =~ /\S\S/) { #diploid site
			$sampleGTcall = &mkDiploidGTcall ($sampleAlleleCalls);
		} else {
			$sampleGTcall = &mkHaploidGTcall ($sampleAlleleCalls);
		}
		if ($sampleGTcall =~ /warning/) {
			print STDERR "$sampleGTcall,$chr-$begin-$end,$idstr,$altstr,$sample,", join ("-", @$sampleAlleleCalls), "\n";
			$sampleGTcall = '.';
			$sampleGTcall = './.' if ($sampleAlleleCalls->[0] =~ /\S\S/);
		}
		$info_NS-- if ($sampleGTcall =~ /\./);
		push (@allGTcalls, $sampleGTcall);
		$info_AF[$1]++ if ($sampleGTcall =~ /^(\d+)/);
		$info_AF[$1]++ if ($sampleGTcall =~ /^\d+\/(\d+)$/);
	}
	$filter = 's50' if ($info_NS/scalar @$samples < 0.5);
	print "$filter\t";
	my ($numAlleles) = 0;
	foreach (@info_AF) {$numAlleles += $_;}
	$numAlleles = 1 if ($numAlleles == 0); #In some loci, all samples are no-called
	my (@alleleFreq) = map (sprintf ('%.2f', $_/$numAlleles), @info_AF);
	shift @alleleFreq; #don't need to print the reference allele frequency
	print "NS=$info_NS;RSID=", join (",", @info_RSid), ";AF=", join (",", @alleleFreq), ";$info_DB\t$format\t", join ("\t", @allGTcalls), "\n";
}

#Convert a list of allele calls to diploid VCF genotype calls
#[11] becomes 1/1
#[01] becomes 0/1
#[00] becomes 0/0
#[11, 00] becomes 1/1
#[01, 01] becomes 1/2
#[00, 01, 1N] becomes 2/3
#[01, 01, 01] becomes ./. (unknown)
#[1N, 00, 00, 00] becomes 1/.
sub mkDiploidGTcall {
	my ($sampleAlleleCalls) = @_;
	#[11, 00, NN, 00], [1N, NN, 01, 00], [1N, NN, 1N, 00], [01, NN, 01, 00], [1N, 00, 00, 00]
	my ($numOfOnes) = 0;
	my ($sampleGTcall) = './.'; #diploid no-call in VCF
	my (%merged);
	for (my $i=0; $i<scalar @$sampleAlleleCalls; $i++) {
		push (@{$merged{$sampleAlleleCalls->[$i]}}, $i+1);
		my ($ones) = $sampleAlleleCalls->[$i] =~ tr/1/1/;
		$numOfOnes += $ones;
	}
	return 'warning-more_than_two_1s' if ($numOfOnes > 2);
	return $merged{'11'}->[0] . '/' . $merged{'11'}->[0] if (exists $merged{'11'});
	if (exists $merged{'01'}) {
		return $merged{'01'}->[0] . '/' . $merged{'01'}->[1] if (defined $merged{'01'}->[1]);
		return $merged{'01'}->[0] . '/' . $merged{'1N'}->[0] if (exists $merged{'1N'});
		return '0/' . $merged{'01'}->[0] if (!exists $merged{'NN'} && !exists $merged{'0N'});
		return $merged{'01'}->[0] . '/.';
	}
	if (exists $merged{'1N'}) {
		return $merged{'1N'}->[0] . '/' . $merged{'1N'}->[1] if (defined $merged{'1N'}->[1]);
		return $merged{'1N'}->[0] . '/.';
	}
	return '0/0' if (!exists $merged{'NN'} && !exists $merged{'0N'});
	return '0/.' if (exists $merged{'0N'} && scalar @{$merged{'0N'}} == 1);
	return $sampleGTcall;
}

#Convert a list of allele calls to haploid VCF genotype calls
#[1] becomes 1
#[N] becomes .
#[0] becomes 0
#[1, 0] becomes 1
#[0, 1, N] becomes 2
#[1, 1] becomes .
sub mkHaploidGTcall {
	my ($sampleAlleleCalls) = @_;
	my ($sampleGTcall) = '.'; #haploid no-call in VCF
	my (%merged);
	for (my $i=0; $i<scalar @$sampleAlleleCalls; $i++) {
		push (@{$merged{$sampleAlleleCalls->[$i]}}, $i+1);
	}
	if (exists $merged{'1'}) {
		return "warning-more_than_one_1_haploid" if (defined $merged{'1'}->[1]);
		return $merged{'1'}->[0];
	}
	return '0' if (!exists $merged{'N'});
	return $sampleGTcall;
}
