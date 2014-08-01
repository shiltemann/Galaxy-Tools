inputfile=system("echo $GPmastervar")
inputfile2=system("echo $GPinputfile2")
xbegin=system("echo $GPstart")
xend=system("echo $GPend")
dotsize=system("echo $GPdotsize")
picwidth=system("echo $GPpicwidth")
picheight=system("echo $GPpicheight")
dotstyle=system("echo $GPdotstyle")
ptitle=system("echo $GPtitle")
yminval=system("echo $yminval")
ymaxval=system("echo $ymaxval")
mydotcolour=system("echo $dotcolour")
mydotcolour2=system("echo $dotcolour2")
drawagrid=system("echo $GPgrid")
logscale=system("echo $GPlogscale")
marker=system("echo $GPmarker")

set term png size picwidth,picheight
set mxtics 4 #minor xtics

if(drawagrid == 1) set grid ytics lc rgb "#bbbbbb" lw 1 lt 0

if(marker) set arrow from marker,graph(0,0) to marker,graph(1,1) nohead

### SNPs ###
unset xlabel
set ytics font "arial,14"
set format x "%.00s Mb" 
set xtics font "arial,14"
set title ptitle
set title font "arial,25"

### mastervar ###
set yrange [yminval:ymaxval]
if(logscale == 2) set logscale y 2
if(logscale == 10) set logscale y	
set xrange [xbegin:xend]
set output "out.png"
plot inputfile u 2:4 with points lc rgb mydotcolour pt dotstyle ps dotsize notitle, inputfile2 u 2:4 with points lc rgb mydotcolour2 pt dotstyle ps dotsize notitle



