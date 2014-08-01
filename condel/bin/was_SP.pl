#!/soft/bin/perl

use strict;

#############################################################################################
# Reading configuration specifics

if (@ARGV < 2){
  die "USAGE: was_SP.pl condel_config_dir methods_score_file\n\n";
}

my $config_dir = $ARGV[0];
my $config_file = "$config_dir/condel_SP.conf";
my $safe_conf = 0;

our %config;

open(CONF, "<$config_file") || die "Could not open $config_file";
my @conf = <CONF>;
for (my $i = 0; $i < @conf; $i++){
  if ($conf[$i] =~ /condel\.dir=\'(\S+)\'/){
    $config{'condel.dir'} = $1;
    $safe_conf++;
  }
  elsif ($conf[$i] =~ /(cutoff\.HumVar\.\w+)=\'(\S+)\'/){
    $config{"$1"} = $2;
    $safe_conf++;
  }
  elsif ($conf[$i] =~ /(max\.HumVar\.\w+)=\'(\S+)\'/){
    $config{"$1"} = $2;
    $safe_conf++;
  }
}

if ($safe_conf < 3){
  die "Malformed config file!!!\n\n";
}

# End configuration specifics
#################################################################################################

#################################################################################################
# Reading complementary cumulative distributions files

my %sift;
my %polyphen;

my %WAS;
my %class;

open(SIFT, $config{'condel.dir'}."/methdist/sift.data");
my @sift = <SIFT>;
close SIFT;

for (my $i = 0; $i < @sift; $i++){
  if ($sift[$i] =~ /(\S+)\s+(\S+)\s+(\S+)/){
    $sift{'tp'}{"$1"} = $2;
    $sift{'tn'}{"$1"} = $3;
  }
}

open(POLYPHEN, $config{'condel.dir'}."/methdist/polyphen.data");
my @polyphen = <POLYPHEN>;
close POLYPHEN;

for (my $i = 0; $i < @polyphen; $i++){
  if ($polyphen[$i] =~ /(\S+)\s+(\S+)\s+(\S+)/){
    $polyphen{'tp'}{"$1"} = $2;
    $polyphen{'tn'}{"$1"} = $3;
  }
}

# Files read
#################################################################################################

my $methscores = $ARGV[1];

#################################################################################################
# Opening and checking input file

open(SCORES, $methscores);
my @score = <SCORES>; 

if ($score[1] !~ /(\S+)\s+(\S+)\s+(\S+)/){
  die "Malformed scores file!!!\n\n";
}

#
##################################################################################################
# Computing WAS

for (my $i = 0; $i < @score; $i++){
  if ($score[$i] =~ /(\S+)\s+(\S+)\s+(\S+)/){
    next if ($1 =~ /id/);
    $WAS{'id'} = $1;
    $WAS{'sift'} = $2;
    $WAS{'pph2'} = $3;
   
    my $base = 0;
    my $int_score = 0;
    my $sift_score = sprintf("%.3f", $WAS{'sift'});
    if ($sift_score <= $config{'cutoff.HumVar.sift'}){
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tn'}{"$sift_score"}));
      $base += 1-$sift{'tn'}{"$sift_score"};
      $class{'sift'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", (1 - $sift_score/$config{'max.HumVar.sift'})*(1-$sift{'tp'}{"$sift_score"}));
      $base += 1-$sift{'tp'}{"$sift_score"};
      $class{'sift'} = 'neutral';
    }
    my $polyphen_score = sprintf("%.3f", $WAS{'pph2'});
    if ($polyphen_score >= $config{'cutoff.HumVar.polyphen'}){
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tn'}{"$polyphen_score"}));
      $base += 1-$polyphen{'tn'}{"$polyphen_score"};
      $class{'polyphen'} = 'deleterious';
    }
    else {
      $int_score += sprintf("%.3f", $polyphen_score/$config{'max.HumVar.polyphen'}*(1-$polyphen{'tp'}{"$polyphen_score"}));
      $base += 1-$polyphen{'tp'}{"$polyphen_score"};
      $class{'polyphen'} = 'neutral';
    }
    if ($base == 0){
      $int_score = -1;
      $class{'condel'} = 'not_computable_was';
    }
    else {
      $int_score = sprintf("%.3f", $int_score/$base);
    }

    if ($int_score >= 0.469){
      $class{'condel'} = 'deleterious';
    }
    elsif ($int_score > 0 && $int_score < 0.469) {
      $class{'condel'} = 'neutral';
    }
    
    print $WAS{'id'}, "\t", $int_score, "\t", $class{'condel'} ,"\n";

  }
}


