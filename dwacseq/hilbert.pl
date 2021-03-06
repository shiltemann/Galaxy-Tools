#!/usr/bin/perl -w

use strict;

use Math::Curve::Hilbert;
#use CGI;
use GD;
use POSIX qw(tmpnam ceil floor);
use File::Copy qw(copy);
use constant LEVEL => 9; # 2^level by 2^level image

my $hilbert_size = 512*512;

#my $query = new CGI;
#my $upload_filehandle = $ARGV[0];
my $inputfile=$ARGV[0];
my $legendfile=$ARGV[1];
my $legendoriginal=$ARGV[2];


copy($legendoriginal,$legendfile);



my $tmp_filename = "file";
my $tmp_picturename = "picture.png";
#print "Content-type:text/html;\n\n", ">/srv/www/genomes.nl/tools/dwac$tmp_filename";
#exit;


print "inputfilee: ";
print $inputfile;
print "\n";
print "tmp_filename: ";
print $tmp_filename;
print "\n";
print "picturename: ";
print $tmp_picturename;
print "\n";
print "legend original: ";
print $legendoriginal;
print "\n";

open FILE, ">", $tmp_filename;


open FILO, "<", $inputfile;

my @data = ();
my @data2 = ();

while ( <FILO> )  
{  
	my @line = split("\t", $_);
	push (@data, $line[3]);
	chomp $line[-1];
	push @data2, join( ' ', @line);
}

# determine data traits

#size
unless ( @data ) {
    print "Content-type:text/html;\n\nNo data found. Please check for using the proper data format!";
    exit;
}
my $nOfDataPoints = scalar @data;

#weight of point
my $low_scale = floor($hilbert_size/$nOfDataPoints);
my $high_scale = ceil($hilbert_size/$nOfDataPoints);

if ($low_scale == 0) {$low_scale = 1}

#median
my @dataSorted = sort { $a <=> $b } @data; 

my $dataMedian;
if ($nOfDataPoints % 2){
	$dataMedian = $dataSorted[int($nOfDataPoints/2)];
}
else{
	$dataMedian = ($dataSorted[$nOfDataPoints/2] + $dataSorted[$nOfDataPoints/2 - 1]) / 2; 
}

print FILE "$dataMedian";
close FILE;

###

my $dim = 2**LEVEL;
my $nOfHilbertPoints = $dim * $dim;

my $hilbert = Math::Curve::Hilbert->new( direction => 'up', max => LEVEL, clockwise => 1, step => 1);

# create a new image and define colors

my $im = new GD::Image($dim+2, $dim+2);

my $black = $im -> colorAllocate(0,0,0);# no data
my $blue = $im -> colorAllocate(0,0,255);# 0 copies
my $skyBlue = $im -> colorAllocate(0,191,255);# 1 copy
my $white = $im -> colorAllocate(255,255,255);# 2 copies
my $pink = $im -> colorAllocate(255, 192, 203);# 3 copies
my $red = $im -> colorAllocate(255,0,0);# 4 copies
my $darkRed = $im -> colorAllocate(139,0,0);# > 4

# determine the number of data points in one pixel
my $pixelResolution = int($nOfDataPoints/$nOfHilbertPoints);
if ($pixelResolution == 0) {$pixelResolution = 1;}

# create an array of bins colors
my @dataCategorized;
if ($nOfHilbertPoints > $nOfDataPoints) {$nOfHilbertPoints = $nOfDataPoints}

for(my $i = 1; $i < $nOfHilbertPoints; $i++) {
	my @slice = ();
	for(my $j = 1; $j < ($pixelResolution+1); $j++) {push(@slice,shift(@data))};
	my $sliceAverage = average(\@slice);
	
	my $binLabel;
	
	if (($sliceAverage >= 0) & ($sliceAverage < $dataMedian/4)) {$binLabel = 1;}
	if (($sliceAverage >= $dataMedian/4) & ($sliceAverage < (3/4)*$dataMedian)) {$binLabel = 2;}
	if (($sliceAverage >= (3/4)*$dataMedian) & ($sliceAverage < (5/4)*$dataMedian)) {$binLabel = 3;}
	if (($sliceAverage >= (5/4)*$dataMedian) & ($sliceAverage < (7/4)*$dataMedian)) {$binLabel = 4;}
	if (($sliceAverage >= (7/4)*$dataMedian) & ($sliceAverage < (9/4)*$dataMedian)) {$binLabel = 5;}
	if ($sliceAverage >= (9/4)*$dataMedian){$binLabel = 6;}

	for (my $counter=1; $counter <= $low_scale; $counter++) { push (@dataCategorized, $binLabel) }
}

# and the last point 
my $sliceAverage = average(\@data);

my $binLabel;
if (($sliceAverage >= 0) & ($sliceAverage < $dataMedian/4)) {$binLabel = 1;}
if (($sliceAverage >= $dataMedian/4) & ($sliceAverage < (3/4)*$dataMedian)) {$binLabel = 2;}
if (($sliceAverage >= (3/4)*$dataMedian) & ($sliceAverage < (5/4)*$dataMedian)) {$binLabel = 3;}
if (($sliceAverage >= (5/4)*$dataMedian) & ($sliceAverage < (7/4)*$dataMedian)) {$binLabel = 4;}
if (($sliceAverage >= (7/4)*$dataMedian) & ($sliceAverage < (9/4)*$dataMedian)) {$binLabel = 5;}
if ($sliceAverage >= (9/4)*$dataMedian){$binLabel = 6;}

for (my $counter=1; $counter <= $low_scale; $counter++) { push (@dataCategorized, $binLabel) }
	
	

# now draw it
my $point = 0;
while (@dataCategorized){
	
	my $colorIndex = shift(@dataCategorized);
	my ($x,$y) = $hilbert->CoordinatesFromPoint($point);
	$point ++;
	
	if ($colorIndex == 1) {$im->setPixel($x,$y,$blue);}
	if ($colorIndex == 2) {$im->setPixel($x,$y,$skyBlue);}
	if ($colorIndex == 3) {$im->setPixel($x,$y,$white);}
	if ($colorIndex == 4) {$im->setPixel($x,$y,$pink);}
	if ($colorIndex == 5) {$im->setPixel($x,$y,$red);}
	if ($colorIndex == 6) {$im->setPixel($x,$y,$darkRed);}
}

# Make sure we are writing to a binary stream
#binmode PICTURE;

# Convert the image to PNG and print it to the file PICTURE
#print $im->png;

# To disable buffering of image content.
#    select(STDOUT);
#    $| = 1;
#    undef $/;

#    print "Content-type: image/png\n\n";
#    print $im->png;
#close PICTURE;


open PICTURE, ">", $tmp_picturename;

# Make sure we are writing to a binary stream
binmode PICTURE;

# Convert the image to PNG and print it to the file PICTURE
print PICTURE $im->png;
close PICTURE;

#print "Content-type: text/html\n\n";
#print <<HTML;

#<html>
#<body>


#<h2>$upload_filehandle</h2>
#<h3> <a href = /dwac/index.html> home </a></h3>

#<h3> Click on a point to see its coordinates</h3>
#<h4> $nOfDataPoints points loaded </h4>

#<a href = "/cgi-bin/hilbert2.cgi?$tmp_filename?$low_scale" target = "clicked"><img src=/dwac$tmp_picturename\.png HSPACE=30 ismap = ismap></a>
#<img src=/img/legenda.png>

#<h3> Data median: $dataMedian </h3>

#<iframe name = "clicked" height = "40" width = "500"> </iframe>

#<body>
#<html>

#HTML


#system('cat hilb9.map');

#open F, '/srv/www/htdocs/slavik/hilbert/hilb9.map';
#while (<F>) {
#    if ( m/\?(\d+)\"/ ) {
#        my $dp = $1;
#        my $alt = $data2[$dp];
#        s/(\?$dp\")/$1 title=\"$alt\"/;
#    }    
#    print;
#}
#close F;

#print "</body>";
#print "</html>";



# average subroutine
sub average { 
@_ == 1 or die ('Sub usage: $average = average(\@array);'); 
my ($array_ref) = @_; 
my $sum; 
my $count = scalar @$array_ref; 
foreach (@$array_ref) { $sum += $_; } 
return $sum / $count; 
} 

