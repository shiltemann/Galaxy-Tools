#!/usr/bin/perl -w
use strict;
use IO::Handle;
use Getopt::Long;
use List::Util qw/ min max sum /;
use Math::Complex;

# Last mod 2012/08/08

# Parameters section
my $reads_per_window = 1000;  # Number of reads in a window
my $bootstraps       = 1000;  # Number of bootstraps (the more the better)
my $window           = 100;   # Size of the window (of windows) for segmentation
my $step             = 9;     # Step of slicing window (the less the better)
my $conf_level       = 0.95;  # Confidence level in change points finding

my $cnv_lower_ratio  = 0.8;   # Do not fine-tune if ratio is between _lower_ and _higher_ 
my $cnv_higher_ratio = 1.2;

my $read_length      = 35;   # Define read length of NGS read

# New things
my $genome_wide      = 0;
my $design_file      = q{};
my $window_is_target = q{};
my $normal_ratio     = q{};
my $output_prefix    = q{};
my $skip_ratios      = q{};
my $test_strand      = q{};
my $ref_strand       = q{};
my $output_format    = 'wig';
my $mapq             = 0;


my $samtools_exe     = ''; # E.g. /usr/local/samtools/0.1.15/, empty if in PATH

my $usage = "---=== DWAC-Seq v. 0.7 ===---
Dynamic Window Approach for CNV calling using Sequencing data

	Main options:
        --test <test_bam>         Sorted and indexed BAM file for test set
        --ref <ref_bam>           Same for reference set
        --chr <chromosome_name>   Chromosome name
        --amb ambiguity (optional, default 10)    Mapping ambiguity (max best mapping positions per read)
        --scale <reads> (optional, default 1000)  Number of reads per window

	Output options:
	--output <prefix>  Prefix for output files
	--output_format <format> Output format: wig(default), vcf
	--skip_ratios <ration_from>,<ratio_to> Do not report segments with ratios between given

	Advanced options:
        --no_fine_tuning (the results are not fine tuned if present)
	--genome_wide_normalization (optional, recommended for cancer samples. Default = chromosome-wide)
	--remove_clonal (removes duplicate reads, mapped to the same chromosomal position)
	--fixed_window <size_bp> use fixed-size windows instead 
	--normal_ratio <ratio>, specify normlaization ratio (test/reference)t
	--read_length <size_bp>, skip all reads that are shorter than specified
	--test_strand [watson|crick], consider only one strand (for Strand-Seq data analysis) 
	--ref_strand [watson|crick], consider only one strand (for Strand-Seq data analysis)
	
	Options for exome enrichment sets
	--design_file <BED_file> (optional, only considers ranges given in the design file)
	--window_is_target (optional)
	Note: design_file option can be combined with fixed/dynamic window approaches
";

my $n_parameters = @ARGV;
die($usage) unless $n_parameters >= 3;
STDOUT->autoflush(1);
warn "Normal_ratio:$normal_ratio\n";

my $no_fine_tuning = 0;  # by default do fine tuning
my $ambiguity      = 0; # default ambiguity
my $remove_clonal  = 0;  # Do not remove clonal by default
my $no_gaps 	   = 0;  # Reads may have gaps when aligned to reference, by default
my $no_second_best = 0;  # Reads may have second best hits when aligned to reference, by default

my ( $test_bam, $ref_bam, $target_chr, $fixed_window );

GetOptions ( 'no_fine_tuning+'            => \$no_fine_tuning,
             'test=s'                     => \$test_bam,
             'ref=s'                      => \$ref_bam,
             'chr=s'                      => \$target_chr,
             'amb:i'                      => \$ambiguity,
             'no_second_best+'            => \$no_second_best,
             'no_gaps+'                   => \$no_gaps,
             'scale:i'                    => \$reads_per_window,
             'mapq+'			  => \$mapq,
             'skip_ratios=s'		  => \$skip_ratios,
             'test_strand=s'		  => \$test_strand,
             'ref_strand=s'		  => \$ref_strand,
             'genome_wide_normalization+' => \$genome_wide,
             'design_file=s'              => \$design_file,
             'window_is_target+'          => \$window_is_target,
             'remove_clonal:i'            => \$remove_clonal,
             'normal_ratio:f'             => \$normal_ratio,
             'fixed_window:i'             => \$fixed_window,
             'read_length:i'              => \$read_length,
             'output=s'                   => \$output_prefix,
             'output_format=s'            => \$output_format,
           );

die("no chromosome specified") unless $target_chr;
if ( $skip_ratios ) {
    die "Incorrect format of skip_ratios tag: $skip_ratios, example: 0.8,1.2" unless $skip_ratios =~ m/^([\.\d]+)\,([\.\d]+)/;
    $cnv_lower_ratio = $1;
    $cnv_higher_ratio = $2;
    die "skip_ratios:Left ratio should not exceed right ratio" if $cnv_lower_ratio > $cnv_higher_ratio;
}
die "Unknown test_strand option $test_strand" if $test_strand and $test_strand !~ m/^(watson|crick)$/;
die "Unknown ref_strand option $ref_strand" if $ref_strand and $ref_strand !~ m/^(watson|crick)$/;

my %design = ();
if ( $design_file ) {
    warn "Reading the design file\n";
    open( F, $design_file ) or die "Cannot read $design_file\n";
    while ( <F> ) {
        chomp;
        my ( $chr, $start, $end ) = split /\t/;
        die "duplicate entries in design file: $chr $start\n" if exists( $design{$chr}{$start});
        next if exists( $design{$chr}{$start});
        $design{$chr}{$start} = $end;
    }
    close F;
    warn "Checking the design file\n";
    foreach my $chr ( keys %design ) {
        my $prev_end = 0;
        foreach my $start ( sort {$a<=>$b} keys %{$design{$chr}} ) {
            die "record $chr\:$start\-$design{$chr}{$start} overlaps with previous segment, ending on $prev_end, please re-cluster\n" if $prev_end >= $start;
            $prev_end = $design{$chr}{$start};
        }
        warn ' ', scalar keys %{$design{$chr}}, " targets on chr $chr\n";
    }
}


# Extract chromosome $target_chr from the BAM files and create output file

my ( $test, $ref ) = ( $test_bam, $ref_bam );
$test =~ s/^.*\///g;
$ref  =~ s/^.*\///g;
$test =~ s/^([a-z\d]+).*$/$1/ if $test =~ m/^[a-z\d]+/;
$ref  =~ s/^([a-z\d]+).*$/$1/ if $ref  =~ m/^[a-z\d]+/;

my $output_name = $output_prefix ? $output_prefix :
                  join( '_', 'tst', $test,
                             'ref', $ref,
                             'amb', $ambiguity,
                             'win', $reads_per_window,
#                             'TEST'
                      );

my @ref   = ();
my @refgc = ();
my $max_coord = 1;
my $last_pos = 0;
my $clonal = 0;
my $clonal_again = 0;
warn "Reading REF amb $ambiguity nsb $no_second_best ng $no_gaps\n";
open( F, $samtools_exe.'samtools view '.$ref_bam.' '. $target_chr.' |') or die "No reference BAM given";
warn $samtools_exe.'samtools view '.$ref_bam.' '. $target_chr.' |'."\n";
my $i = 0;
my @starts = keys %design ? sort {$a<=>$b} keys %{$design{$target_chr}} : ();
while ( <F> ) {
    last if $design_file and !@starts;
    next if $ambiguity and m/X0\:i\:(\d+)/ and $1 > $ambiguity;
    next if $no_second_best and m/X1\:i\:(\d+)/ and $1 > 0;
    next if $no_gaps and m/XO\:i\:(\d+)/ and $1 > 0;
    next if $no_gaps and m/XG\:i\:(\d+)/ and $1 > 0;
    
    my @values = split /\t/;
    next if $mapq and $values[4]< $mapq;
    next if $remove_clonal and $values[1]&1024;
    next if length($values[9]) < $read_length;
    if ( $remove_clonal and $values[3] == $last_pos ) {
        $clonal++;
        next if $remove_clonal > 0 and ++$clonal_again > $remove_clonal;
    } elsif ( $remove_clonal and $values[3] != $last_pos ) {
        $clonal_again = 0;
    }
    next if $ref_strand eq 'watson' and ( ( $values[1] & 16 and $values[1] & 64 ) or ( not($values[1] & 16 ) and $values[1] & 128 ) );
    next if $ref_strand eq 'crick'  and ( ( $values[1] & 16 and $values[1] & 128 ) or ( not($values[1] & 16) and $values[1] & 64  ) );
    $max_coord = $values[3];
    my $gc = $values[9];
    #warn "gc $gc\n";
    $gc =~ s/[^GATCgtac]//ig;
    my $len = length( $gc );
    next unless $len;
    $gc=~s/[^GCgc]//ig;
    if ( $design_file ) {
        while ( $design{$target_chr}{$starts[0]} < $values[3] ) {
            shift @starts;
            last unless @starts;
        }
        if ( @starts and $starts[0] <= $values[3] and $values[3] <= $design{$target_chr}{$starts[0]} ) {
            push @ref, $values[3];
            push @refgc, length($gc)/$len;
        } 
    }
    else {
        push @ref, $values[3];
        push @refgc, length($gc)/$len;
    }
    $last_pos = $values[3];
    print "\r",$i/1_000_000,"M reads " if ++$i % 100_000 == 0;
}
print "Total: $i raw ", scalar @ref, " target reads                \n";
print "Filtered out  $clonal clonal; reads                         \n" if $remove_clonal;
close F;

die("No reference reads on this chromosome or not sorted BAM file") unless $#ref > 0;
my @test = ();
warn "Reading TEST\n";
$i = 0;
$last_pos = 0;
$clonal = 0;
@starts = keys %design ? sort {$a<=>$b} keys %{$design{$target_chr}} : ();
open( F, $samtools_exe.'samtools view '.$test_bam.' '. $target_chr.' |') or die "No test BAM given";
warn $samtools_exe.'samtools view '.$test_bam.' '. $target_chr.' |'."\n";
while ( <F> ){
    last if $design_file and !@starts;
    next if $ambiguity and m/X0\:i\:(\d+)/ and $1 > $ambiguity;
    next if $no_second_best and m/X1\:i\:(\d+)/ and $1 > 0;
    next if $no_gaps and m/XO\:i\:(\d+)/ and $1 > 0;
    next if $no_gaps and m/XG\:i\:(\d+)/ and $1 > 0;
    my @values = split /\t/;
    next if $mapq and $values[4]< $mapq;
    next if $remove_clonal and $values[1]&1024;
    next if length($values[9]) < $read_length;
    if ( $remove_clonal and $values[3] == $last_pos ) {
        $clonal++;
        next if $remove_clonal > 0 and ++$clonal_again > $remove_clonal;
    } elsif ( $remove_clonal and $values[3] != $last_pos ) {
        $clonal_again = 0;
    }
    next if $test_strand eq 'watson' and ( ( $values[1] & 16 and $values[1] & 64 ) or ( not($values[1] & 16 ) and $values[1] & 128 ) );
    next if $test_strand eq 'crick'  and ( ( $values[1] & 16 and $values[1] & 128 ) or ( not($values[1] & 16) and $values[1] & 64  ) );
    if ( $design_file ) {
        while ( $design{$target_chr}{$starts[0]} < $values[3] ) {
            shift @starts;
            last unless @starts;
        }
        push @test, $values[3] if @starts and $starts[0] <= $values[3] and $values[3] <= $design{$target_chr}{$starts[0]};
    }
    else {
        push @test, $values[3];
    }
    $last_pos = $values[3];
    print "\r",$i/1_000_000,"M reads " if ++$i % 100_000 == 0;
}
print "Total: $i raw ",scalar @test, " target reads                \n";
print "Filtered out  $clonal clonal; reads                         \n" if $remove_clonal;
close F;

die("No test reads on this chromosome or not sorted BAM file") unless $#test > 0;

# step 1 - counting number of reads per window

my $window_start = $ref[0];
my $window_end   = -1;
my $read;
my $test_reads = 0;
my $ref_reads  = 0;
#


# Setting windows' boundaries
my %windows = ();
# Mod June 1st nado proveryat# consistentnost# opciy i delat# merge design listu
if ( $window_is_target ) {
    warn "Setting boundaries of windows, window is target\n";
    foreach my $start ( sort {$a<=>$b} keys %{$design{$target_chr}} ) {
        $windows{$start} = $design{$target_chr}{$start};
    }
}
elsif ( $fixed_window and not(keys %design) ) {
    warn "Setting boundaries of windows, window is $fixed_window bp\n";
    my $start = 1;
    while ( $start + $fixed_window < $max_coord ) {
        $windows{$start} = $start + $fixed_window - 1;
        $start+=$fixed_window; 
    }
}
elsif ( $fixed_window and keys %design ) { # also possible
    warn "Setting boundaries of windows, window is $fixed_window bp across target regions\n";
    my $accu = 0;
    my $window_start;
    foreach my $start ( sort {$a<=>$b} keys %{$design{$target_chr}} ) {
        warn "Started with $start\n";
        $window_start = $start if $accu == 0;
        my $end = $design{$target_chr}{$start};
        warn "Started with target $start $end\n";
        my $window_end;
        my $current_start = $start;
        while ( $accu + $end - $current_start >= $fixed_window ) {
            $window_end = $current_start + $fixed_window - $accu - 1;
            $windows{$window_start} = $window_end;
            warn "Placed window $window_start $window_end\n";
            $window_start = $window_end + 1;
            $current_start = $window_start;
            $accu = 0;
        }
        if ( $accu == 0 ) {
            $accu = $end - $window_end;
        }
        else {
            $accu += ( $end - $start + 1 );
        }
    }
}
else {
    warn "Setting boundaries of windows, window is dynamic size, $reads_per_window reads per window\n";
    foreach my $ref_pos ( @ref ) {
        die "Reference BAM is not sorted" if ($ref_pos + 1) < $window_start; 
        $ref_reads++ if $ref_pos >= $window_start;
        if ( $ref_reads >= $reads_per_window ) {
            $window_end = $ref_pos;
            $windows{$window_start} = $window_end;
            $window_start = $ref_pos + 1;
            $ref_reads = 0;
        }
    }
}

warn "Counting reads in windows for test set\n";
my %test_win  = ();
@starts = sort {$a<=>$b} keys %windows;
foreach my $test_pos ( @test ) {
    last unless @starts;
    if ( $test_pos < $starts[0] ) {
        ;
    } 
    elsif ( $test_pos <= $windows{$starts[0]} ) {
        $test_win{$starts[0]}++;
    }
    else {
        shift @starts;
        redo;
    }
}


#Counting number of reads in a reference window
warn "Counting reads in windows for reference set\n";
my %ref_win = ();
my @full = ();
@starts = sort {$a<=>$b} keys %windows;
foreach my $ref_pos ( @ref ) {
    next unless @starts;
    if ( $ref_pos < $starts[0] ) {
        ;
    } 
    elsif ( $ref_pos <= $windows{$starts[0]} ) {
        $ref_win{$starts[0]}++;
    }
    else {
        $test_win{$starts[0]} = 0 unless exists($test_win{$starts[0]});
        unless ( exists($ref_win{$starts[0]}) ) {
            delete $windows{$starts[0]};
            delete $test_win{$starts[0]};
        }
        else {
            push @full, $test_win{$starts[0]} / $ref_win{$starts[0]};
        }
        shift @starts;
        redo;
    }
}



my %ref_gc_win = ();
@starts = sort {$a<=>$b} keys %windows;
foreach my $ele ( 0 .. $#ref ) {
    last unless @starts;
    my $pos = $ref[$ele];
    if ( $pos <= $windows{$starts[0]} ) {
        push @{$ref_gc_win{$starts[0]}}, $refgc[$ele];
    }
    else {
        shift @starts;
        redo;
    }
}

my @gc_norm = ();
my @ratio_norm = ();
my %gc_win = ();
foreach my $start ( keys %ref_gc_win ) {
    my ($gc, $i) = (0,0);
    map { $gc+=$_ ; $i++} @{$ref_gc_win{$start}};
    $gc /= $i;
    push @gc_norm, $gc;
    $gc_win{$start} = $gc;
    push @ratio_norm, $test_win{$start}/$ref_win{$start};
}

open( F, '>', $output_name.'_chr'.$target_chr.'.nrm') or die "Cannot write on disk";
my @gc_buf = ();
my @ratio_buf = ();
my %gc2ratio = ();
$i = 0;
foreach my $ele ( sort { $gc_norm[$a]<=>$gc_norm[$b] } (0..$#gc_norm) ) {
    push @gc_buf, $gc_norm[$ele];
    push @ratio_buf, $ratio_norm[$ele];
    shift @gc_buf if @gc_buf > 1000;
    shift @ratio_buf if @ratio_buf > 1000;
    if ( ++$i % 100 == 0 and @gc_buf == 1000 ) {
        my $ratio = 0;
        map { $ratio += $_ } @ratio_buf;
        print F $gc_buf[499], "\t", $ratio/1000, "\n";
        $gc2ratio{$gc_buf[499]} = $ratio/1000;
    } 
}

# Normalize for GC
#foreach my $start ( keys %test_win ) {
#    my $gc = $gc_win{$start};
#    my $best_gc = 100;
#    foreach my $gc2 ( keys %gc2ratio ) {
#         $best_gc = $gc2 if abs( $gc2 - $gc ) < abs( $best_gc - $gc );
#    }
#    $test_win{$start} = int( $test_win{$start} /$gc2ratio{$best_gc} + 0.5 );
#} 



my $n_windows = scalar keys %windows;
@starts = sort {$a<=>$b} keys %windows;
open( F, '>', $output_name.'_chr'.$target_chr.'.win') or die "Cannot write on disk";
foreach my $i( 0 .. $n_windows-1 ) {
     print F  join( "\t", $target_chr, $starts[$i], $windows{$starts[$i]}, $test_win{$starts[$i]}, $ref_win{$starts[$i]} ),"\n";
}
close F;


# Step 2 - segmentation
warn "Segmenting coverage profile\n";
open( F, '>', $output_name.'_chr'.$target_chr.'.brk') or die "Cannot write on disk";
my $data_size = scalar keys %windows;
@starts = sort {$a<=>$b} keys %windows;
my %change_points = ();
my $progress_prev = 0;

for ( my $counter = 0; $counter < $data_size - 1; $counter += $step ) {
    print "\rwindow ",$counter+1,' / ', $data_size, ' ';
    my $s = $counter;
    my $e = $counter + $window - 1;
    next if $e > $#starts;
    my ( $point, $conf ) = find_cp( $s, $e );
    if ( $conf >= $conf_level and ( !exists($change_points{ $point }) or $change_points{ $point } < $conf ) ) { # change point is located after window with number $point
        $change_points{ $point } = $conf;
        my $start = $starts[$point];
#        print "Found breakpoint at $windows{$start} confidence $conf\n";
        print F join( "\t", $target_chr, $windows{$start}, $conf ), "\n";
    }
}
close F;

## add last point
$change_points{ $data_size-1} = 1;

# Averaging of copy-number for each segment
print "\nAveraging ratio between breakpoints\n";
my @reads_per_window_av = ();
my $previous_change_point = -1; # before window 0
foreach  my $current_change_point ( sort { $a <=> $b } keys %change_points ) {
    # average from previous change point
    my @values = ();
    foreach my $k ( ($previous_change_point + 1) .. $current_change_point ) {
        die "Attempt to get over end of chromosome" if $k > $#starts;
        $test_win{$starts[$k]} = 0 unless exists($test_win{$starts[$k]});
        unless ( $ref_win{$starts[$k]} ) { warn "Zero ref reads: $starts[$k]\n"; next }
        push @values, $test_win{$starts[$k]} / $ref_win{$starts[$k]};
    }
    my $median = median(\@values);
    # save averaged data
    for my $k ( ( $previous_change_point + 1) .. $current_change_point ) {
        $reads_per_window_av[$k] = $median;
    }
    $previous_change_point = $current_change_point;
}

my $med = $normal_ratio ? $normal_ratio :
          $genome_wide  ? median_whole_genome() :
          median(\@reads_per_window_av); 
warn "Using median ratio: $med\n";

print "Correcting ratios using median ratio ($med)\n";
my @reads_normalized = ();
foreach ( @reads_per_window_av ) { 
    my $rounded =  sprintf ( "%.3f", $_/$med );
    push ( @reads_normalized, $rounded );
}

print "Merging adjacent regions that have similar copy number\n";
my %data_out = ();
my $r = $reads_normalized[0];
my $be = 0;
my $e;
foreach  my $counter ( 1 .. $n_windows - 1 ) {
    next unless abs( $reads_normalized[$counter] - $r ) > 0.1;
    $e = $counter;
    $data_out{ $be."\t".($e - 1) } = $r > 0 ? $r : 0.00001; #otherwise logn will fail
    $r = $reads_normalized[$counter];
    $be = $counter;
}
# ... and the last segment
$data_out{ $be."\t". $#starts} = $r > 0 ? $r : 0.00001;

print "Saving results\n";

my %data_out_sorted;

foreach (keys %data_out) {
    my ( $s, $e ) = split /\t/;
    $data_out_sorted{$s} = join("\t", $data_out{ $s."\t".$e}, $starts[$s], $windows{$starts[$e]});
}

open( F, '>', $output_name.'_chr'.$target_chr.'.seg') or die "Cannot write on disk";
foreach (sort{$a<=>$b} keys %data_out_sorted){
    print F "$target_chr\t$data_out_sorted{$_}\n";
}
close F;

open( F, '>', $output_name.'_chr'.$target_chr.'.avg') or die "Cannot write on disk";
foreach my $i( 0 .. $n_windows-1 ) {
     print F join( "\t", $target_chr, $starts[$i], $windows{$starts[$i]}, $reads_per_window_av[$i]),"\n";
}
close F;

# Do fine-tuning when requested
exit if $no_fine_tuning;
warn "Fine-tuning CNV boundaries\n";

my %seen = ();

## process windows sizes 
open( F, '>', $output_name.'_chr'.$target_chr.'.txt') or die "Cannot write on disk";
open( FD, '>', $output_name.'_chr'.$target_chr.'.diag') or die "Cannot write on disk";
my $nOfSegments = scalar keys %data_out;
foreach my $range ( sort { abs(logn($data_out{$b},2)) <=> abs(logn($data_out{$a},2)) } keys %data_out ) {
    my ( $swin, $ewin ) = split /\t/, $range;
    my $r = $data_out{$range}; # Starting ratio 
    next unless $r <= $cnv_lower_ratio or $r >= $cnv_higher_ratio; 

    # Search space for optimal start 
    my $s1 = $swin > 0 ? $starts[ $swin-1 ]: $starts[0];
    my $s2 = $starts[ $swin ];
    my $s3 = $windows{ $starts[$swin] };

    # Search space for optimal end
    my $e1 = $starts[ $ewin ];
    my $e2 = $ewin < $#starts ? $starts[ $ewin + 1 ] : $windows{ $starts[ $ewin ] };
    my $e3 = $ewin < $#starts ? $windows{ $starts[ $ewin + 1 ] } : $windows{ $starts[ $ewin ] };
    
    # Local normalization coeff
    my ( $local_gc, $local_cnt ) = ( 0, 0 );
    foreach my $start ( keys %gc_win ) {
        next if $start < $s2 or $start > $e2;
        $local_gc += $gc_win{$start};
        $local_cnt++;
    }
    $local_gc /= $local_cnt;
    my $best_gc = 100;
    foreach my $gc2 ( keys %gc2ratio ) {
         $best_gc = $gc2 if abs( $gc2 - $local_gc ) < abs( $best_gc - $local_gc );
    }
    my $local_norm = $gc2ratio{$best_gc};

    print FD join( "\t", 'before', $s1, $s2, $s3, 'to', $e1, $e2, $e3, 'ratio', $r), "\n";
    print FD join( "\t", 'reads count test, ref:', calculate_covs( $s2, $e2 ),  'local_norm:', $local_norm ), "\n";
    print "Fine tuning $s2 ($s1\.\.$s3)  $e2 ($e1\.\.$e3) ratio $r\n";


#    print "Fine-tuning $s1 / $s3 / $e1 with $s2 / $e3 / $e2 (was $r)\n";
    my %from_s1 = ();

    print "\rPreparing finetuning for ref set           ";
    my %ref_from_s1 = (); 
    my $clone_count = 0;
    foreach my $pos ( @ref ) {
        next if $pos < $s1;
        last if $pos > $e3;
        $clone_count++;
        $ref_from_s1{$e1} = $clone_count if $e1 >= $pos;
        next unless ( $s1 <= $pos and $pos <= $s3 ) or ( $e1 <= $pos and $pos <= $e3 );
        $ref_from_s1{$pos} = $clone_count;
        $from_s1{$pos} = 1;
    }

    my %test_from_s1 = ();
    print "\rPreparing finetuning for test set           ";
    $clone_count = 0;
    foreach my $pos ( @test ) {
        next if $pos < $s1;
        last if $pos > $e3;
        $clone_count++;
        $test_from_s1{$e1} = $clone_count if $e1 >= $pos;
        next unless ( $s1 <= $pos and $pos <= $s3 ) or ( $e1 <= $pos and $pos <= $e3 );
        $test_from_s1{$pos} = $clone_count;
        $from_s1{$pos} = 1;
    }
    $from_s1{$e1} = 1;

    my ( $ref_clones, $test_clones ) = ( 0, 0 );
    print "\rMerging ref and test              ";
    foreach my $pos ( sort {$a<=>$b} keys %from_s1 ) {
        if ( exists( $ref_from_s1{$pos} ) ) {
            $ref_clones = $ref_from_s1{$pos}
        }
        else {
            $ref_from_s1{$pos} = $ref_clones;
        }
        if ( exists( $test_from_s1{$pos} ) ) {
            $test_clones = $test_from_s1{$pos}
        }
        else {
            $test_from_s1{$pos} = $test_clones;
        }
    }

    print "\rFinding best range                                     ";
    my ( $best_score, $best_s, $best_e, $best_size, $best_r, $best_test, $best_ref ) = ( -10_000, $s2, $e2, $e2-$s2+1, $r, 0, 0 );

    my $test_count = '?';
    my $ref_count  = '?';
    my $init_score = $r;
    
    my @ends = sort {$a<=>$b} grep { $e1 <= $_ and $_ <= $e3 } keys %from_s1;
    foreach my $s ( sort {$a<=>$b} grep { $s1 <= $_ and $_ <= $s3 } keys %from_s1 ) {
SEGEND:
        foreach my $e ( @ends ) {
            next if $s >= $e;
            foreach my $pos ( keys %seen ) {
                next SEGEND unless $pos >= $e or $s >= $seen{$pos};
            }
            die "Found no test from e$e" unless exists($test_from_s1{$e} );
            die "Found no test from s$s" unless exists($test_from_s1{$s} );
            die "Found no ref from e$e" unless exists($ref_from_s1{$e} );
            die "Found no ref from s$s" unless exists($ref_from_s1{$s} );

            $test_count = $test_from_s1{$e} - $test_from_s1{$s};
            $ref_count  = $ref_from_s1{$e} - $ref_from_s1{$s};

            $test_count = 0.1 if $test_count==0;
            $ref_count  = 0.1 if $ref_count==0;
            # Avoid sub-window calls
            next if $r > 1 and $test_count < $reads_per_window;
            next if $r < 1 and $ref_count < $reads_per_window;


            my $score = $r > 1 ? ( $test_count/10000 + $test_count / $ref_count ) : ( $ref_count/10000 + $ref_count / $test_count );
            if ( !$best_score or $best_score < $score or ( $best_score == $score and $best_size < ( $e - $s +1 ) ) ) {
                ( $best_score, $best_s, $best_e, $best_size, $best_r, $best_test, $best_ref ) = ( $score, $s, $e, $e - $s + 1, $test_count / $ref_count / $local_norm, $test_count, $ref_count ) ;
            }
        } 
    }
    warn "test $test_count ref $ref_count best $best_r norm $local_norm\n";
    
    if ( $best_score < 0 or $ref_count < $reads_per_window ) {
        print FD "Suboptimal call, Best score: $best_score ref_count: $ref_count\n";
        next;
    }

    # moved here
#    my $overlaps = 0;
#    foreach my $pos ( keys %seen ) {
#        $overlaps = 1 unless $pos > $best_e+$read_length-1  or $best_s > $seen{$pos};
#    }
#    next if $overlaps;

    $best_r = sprintf ( "%.3f", ($best_r) );
    next if $best_r > $cnv_lower_ratio and $best_r < $cnv_higher_ratio;
#    $test_count = $test_from_s1{$best_e} - $test_from_s1{$best_s};
#    $ref_count  = $ref_from_s1{$best_e} - $ref_from_s1{$best_s};
    $seen{$best_s} = $best_e; # to prevent further overlaps

    print FD join( "\t", 'after', $best_s, 'to', $best_e, 'ratio', $best_r), "\n";
    print FD join( "\t", 'reads count test, ref:', calculate_covs( $best_s, $best_e ),  'norm ', $local_norm ), "\n";

#    print join( "\t", "\nTuned: ", $best_s, $best_e+$read_length-1, $best_r, $test_count, $ref_count ), "\n";
    if ( uc($output_format) eq 'VCF' and $best_r >= 1 ) {
        print F join( "\t", $target_chr, $best_s, '<DUP>', 10, 'PASS', 
                            'IMPRECISE;SVTYPE=DUP;END='.($best_e+$read_length-1).
                            ';SVLEN='.($best_e+$read_length-$best_s),
                            'GT:GQ:CN',
                            './.:0:'.$best_r*2,
                     ), "\n";
    }
    elsif ( uc($output_format) eq 'VCF' and $best_r < 1 ) {
        print F join( "\t", $target_chr, $best_s, '<DEL>', 10, 'PASS', 
                            'IMPRECISE;SVTYPE=DEL;END='.($best_e+$read_length-1).
                            ';SVLEN='.($best_e+$read_length-$best_s),
                            'GT:GQ:CN',
                            './.:0:'.$best_r*2,
                     ), "\n";
    }
    else { 
        print F join( "\t", $target_chr, $best_s, $best_e+$read_length-1, $best_r, $best_test, $best_ref,
#                            $test_from_s1{$best_s}, $test_from_s1{$best_e},
#                            $ref_from_s1{$best_s}, $ref_from_s1{$best_e},
         ), "\n";
    }
}
print F "#END\n";
close F;


### SUBROUTINES
###


# Change Point candidate
sub find_cp {
    my $s1 = shift;
    my $e1 = shift; 
    
    my $size = $e1 - $s1 + 1; 
    return if $size < 3;
    my $m = 0; 
    map { $m += $_ } @full[$s1..$e1];
    $m /= $size; # m = average ratio

    # compute cumulative sums
    my @S = ( $full[$s1] - $m ) ;
    foreach my $counter ( $s1+1 .. $e1 ) {
        push @S, $S[-1] + $full[$counter] - $m;
    }

    map {$_ = abs($_)} @S;
    my $maximum = max(@S);
    my $Sd = $maximum - min(@S);

    my $point = 0; # max point
    while ( $S[$point] < $maximum ) { $point++ } 
    die "Did not reach max" if $S[$point] != $maximum;

    # Bootstrap
    my @B = ();
    foreach  my $strap ( 1 .. $bootstraps ) {
        my @arr = @full[$s1..$e1];
        foreach my $i (reverse ( 1 .. $#arr ) ) {
            my $j = int( rand( $i + 1 ) );
            @arr[$i,$j] = @arr[$j,$i];
        }
        # compute cumulative sums
        @S = ( $arr[0] - $m );
        for my $counter ( 1 .. $#arr ) {
            push @S, $S[-1] + $arr[$counter] - $m;
        }
        map {$_ = abs($_)} @S;
        push @B, (max(@S) - min(@S));
    }

    my $x = 0;
    foreach (@B) { if ($_ < $Sd) {$x++} }

    my $conf = $x/$bootstraps;
    return ($s1 + $point, $conf);
}
# end of Change Point 

sub median { 
    @_ == 1 or die ('Sub usage: $median = median(\@array);'); 
    my ($array_ref) = @_; 
    my $count = scalar @$array_ref; 
    # Sort a COPY of the array, leaving the original untouched 
    my @array = sort { $a <=> $b } @$array_ref; 
    if ( $count % 2 ) { # odd 
        return $array[int($count/2)]; 
    }
    else { # even = avg of 2 middle elements
        return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
    } 
} 

sub median_whole_genome {

    my $median = 1;
    warn "Computing genome-wide median ratio\n";

    # Get the set of all chrmosome names
    my @genome_wide_ratios = (); # Set of ratios over all windows (subject to farther median)
    open( F, $samtools_exe."samtools idxstats $ref_bam |") or die "please update samtools";
    while ( <F> ) {
        next if m/\*/;
        my ( $chr ) = split /\t/;
        print "\rDetermining ratio for chr $chr step 1     ";

        my %gw_windows = ();
        my $prev_start = 0;
        my $last_pos   = 0;
        my $count = 0;
        my @starts = keys %design ? sort {$a<=>$b} keys %{$design{$chr}} : ();
        print "\rDetermining ratio for chr $chr step 1 (define win)     ";
        open( F1, $samtools_exe.'samtools view '.$ref_bam.' '. $target_chr.' |') or die "No reference BAM given";
        while ( <F1> ) {
            last if $design_file and !@starts;
            next unless m/X0\:i\:(\d+)/;
            next unless $1 <= $ambiguity;
            next if $no_second_best and m/X1\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XO\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XG\:i\:(\d+)/ and $1 > 0;
            my @values = split /\t/;
            next if length($values[9]) < $read_length;
            next if $remove_clonal and $values[3] == $last_pos;
            next if $ref_strand eq 'watson' and ( ( $values[1] & 16 and $values[1] & 64 ) or ( not($values[1] & 16 ) and $values[1] & 128 ) );
            next if $ref_strand eq 'crick'  and ( ( $values[1] & 16 and $values[1] & 128 ) or ( not($values[1] & 16) and $values[1] & 64  ) );
            if ( $design_file ) {
                while ( $design{$chr}{$starts[0]} < $values[3] ) {
                    shift @starts;
                    last unless @starts;
                }
                next unless @starts and $starts[0] <= $values[3] and $values[3] <= $design{$chr}{$starts[0]};
            }
            $last_pos = $values[3];
            next if ++$count < $reads_per_window;
            $gw_windows{$prev_start} = $values[3];
            $prev_start = $values[3]+1;
            $count = 0;
        }
        close F1;

        print "\rDetermining ratio for chr $chr step 2 (count test)    ";
        my @gw_starts = sort {$a<=>$b} keys %gw_windows;
        my %gw_count_test = ();
        $last_pos = 0;
        @starts = keys %design ? sort {$a<=>$b} keys %{$design{$chr}} : ();
        open( F1, $samtools_exe.'samtools view '.$test_bam.' '. $target_chr.' |') or die "No reference BAM given";
        while ( <F1> ) {
            last if $design_file and !@starts;
            next unless m/X0\:i\:(\d+)/;
            next if $no_second_best and m/X1\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XO\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XG\:i\:(\d+)/ and $1 > 0;
            next unless $1 <= $ambiguity;
            next unless @gw_starts;
            my @values = split /\t/;
            next if length($values[9]) < $read_length;
            next if $remove_clonal and $values[3] == $last_pos;
            next if $test_strand eq 'watson' and ( ( $values[1] & 16 and $values[1] & 64 ) or ( not($values[1] & 16 ) and $values[1] & 128 ) );
            next if $test_strand eq 'crick'  and ( ( $values[1] & 16 and $values[1] & 128 ) or ( not($values[1] & 16) and $values[1] & 64  ) );
            if ( $design_file ) {
                while ( $design{$chr}{$starts[0]} < $values[3] ) {
                    shift @starts;
                    last unless @starts;
                }
                next unless @starts and $starts[0] <= $values[3] and $values[3] <= $design{$chr}{$starts[0]};
            }
#            next unless @starts;
#            warn "v$values[3] s$starts[0] g$gw_starts[0] w$gw_windows{$gw_starts[0]}\n";
            if ( @starts and $values[3] < $starts[0] ) {
                ;
            } 
            elsif ( $values[3] <= $gw_windows{$gw_starts[0]} ) {
                $gw_count_test{$gw_starts[0]}++;
            }
            else {
                shift @gw_starts;
                redo;
            }
            $last_pos = $values[3];
        }
        close F1;

        print "\rDetermining ratio for chr $chr step 3 (count ref)    ";
        @gw_starts = sort {$a<=>$b} keys %gw_windows;
        my %gw_count_ref = ();
        $last_pos = 0;
        @starts = keys %design ? sort {$a<=>$b} keys %{$design{$chr}} : ();
        open( F1, $samtools_exe.'samtools view '.$ref_bam.' '. $target_chr.' |') or die "No reference BAM given";
        while ( <F1> ) {
            last if $design_file and !@starts;
            next unless m/X0\:i\:(\d+)/;
            next if $no_second_best and m/X1\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XO\:i\:(\d+)/ and $1 > 0;
            next if $no_gaps and m/XG\:i\:(\d+)/ and $1 > 0;
            next unless $1 <= $ambiguity;
            next unless @gw_starts;
            my @values = split /\t/;
            next if length($values[9]) < $read_length;
            next if $remove_clonal and $values[3] == $last_pos;
            next if $ref_strand eq 'watson' and ( ( $values[1] & 16 and $values[1] & 64 ) or ( not($values[1] & 16 ) and $values[1] & 128 ) );
            next if $ref_strand eq 'crick'  and ( ( $values[1] & 16 and $values[1] & 128 ) or ( not($values[1] & 16) and $values[1] & 64  ) );
            if ( $design_file ) {
                while ( $design{$chr}{$starts[0]} < $values[3] ) {
                    shift @starts;
                    last unless @starts;
                }
                next unless @starts and $starts[0] <= $values[3] and $values[3] <= $design{$chr}{$starts[0]};
            }
            if ( $values[3] < $gw_starts[0] ) {
                ;
            } 
            elsif ( $values[3]<= $gw_windows{$gw_starts[0]} ) {
                $gw_count_ref{$gw_starts[0]}++;
            }
            else {
                $gw_count_test{$gw_starts[0]} = 0 unless exists($gw_count_test{$gw_starts[0]}); 
#                warn $gw_count_test{$gw_starts[0]}, ' ', $gw_count_ref{$gw_starts[0]}, "\n";
                push @genome_wide_ratios, $gw_count_test{$gw_starts[0]} / $gw_count_ref{$gw_starts[0]};
                shift @gw_starts;
                redo;
            }
            $last_pos = $values[3];
        }
    } # while reading chromosomes from samtools idxstat   
    my $gw_median = median( \@genome_wide_ratios );
    warn "Genome-wide median :", $gw_median, ' entries: ',scalar(@genome_wide_ratios),"\n";
    return $gw_median;
}

sub calculate_covs {
    my $s = shift;
    my $e = shift;
    return ( scalar( grep {$s <= $_ and $_ <= $e } @test ), scalar( grep {$s <= $_ and $_ <= $e } @ref ) );
}
