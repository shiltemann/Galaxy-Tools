#!/bin/bash

#usage: $0 <mastervarfile> <picture width in pixels> <chromosome> [ <start> <end> ]

echo "command: $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} "
if [[ $# -ne 21 ]]
then
	echo "unexpected number of variables ($#), expected 9. exiting."
	echo "command: $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16}"
	#exit
fi

inputfile=$1
picwidth=$2
picheight=$3
dotsize=$4
chrom=$5
start=$6
end=$7
build=$8
subpages=$9 # Y or N
yminv=${10}
ymaxv=${11}
chrprefix=${12}
dotcol="#${13:6}"
inputfile2=${14}
dotcol2="#${15:6}"
stepsize=${16}
drawgrid=${17}
logscale=${18}
movavg=${19}
marker=${20}
scriptdir=${21}

#echo "dotcol: $dotcol"
if [ $chrom == "chrX" ]
then
	chrom="chr23"
fi
if [ $chrom == "chrY" ]
then
	chrom="chr24"
fi
if [ $chrom == "chrM" ]
then
	chrom="chr25"
fi

if [ $build=="hg18" ]
then
	chrlen[1]=247249719  
	chrlen[2]=242951149  
	chrlen[3]=199501827  
	chrlen[4]=191273063  
	chrlen[5]=180857866  
	chrlen[6]=170899992  
	chrlen[7]=158821424  
	chrlen[8]=146274826  
	chrlen[9]=140273252  
	chrlen[10]=135374737  
	chrlen[11]=134452384  
	chrlen[12]=132349534  
	chrlen[13]=114142980  
	chrlen[14]=106368585  
	chrlen[15]=100338915  
	chrlen[16]=88827254  
	chrlen[17]=78774742  
	chrlen[18]=76117153  
	chrlen[19]=63811651  
	chrlen[20]=62435964  
	chrlen[21]=46944323  
	chrlen[22]=49691432  
	chrlen[23]=154913754  #X
	chrlen[24]=57772954   #Y
	chrlen[25]=16571      #M
fi

if [ $build=="hg19" ]
then
	chrlen[1]=249250621  
	chrlen[2]=243199373  
	chrlen[3]=198022430  
	chrlen[4]=191154276  
	chrlen[5]=180915260  
	chrlen[6]=171115067  
	chrlen[7]=159138663  
	chrlen[8]=146364022  
	chrlen[9]=141213431  
	chrlen[10]=135534747  
	chrlen[11]=135006516  
	chrlen[12]=133851895  
	chrlen[13]=115169878  
	chrlen[14]=107349540  
	chrlen[15]=102531392  
	chrlen[16]=90354753  
	chrlen[17]=81195210  
	chrlen[18]=78077248  
	chrlen[19]=59128983  
	chrlen[20]=63025520  
	chrlen[21]=48129895  
	chrlen[22]=51304566  
	chrlen[23]=155270560 #X 
	chrlen[24]=59373566  #Y 
	chrlen[25]=16571     #M
fi



if [ $start -eq -1 ]
then
 	start=0
fi

if [ $end -eq -1 ]
then
 	m=${chrom:3}	
 	end=${chrlen[$m]} 		
fi

if [ $chrom == "chr23" ]
then
	chrom="chrX"
	m="X"
fi
if [ $chrom == "chr24" ]
then
	chrom="chrY"
	m="Y"
fi
if [ $chrom == "chr25" ]
then
	chrom="chrM"
	m="M"
fi


#get allelefrequencies
#$3 = chromosome
#$4 = begin
#$7 = vartype
#$22 = allele1RC
#$23 = allele2RC
#$25 = totalRC 


export GPmastervar=gnuplot.data
export GPinputfile2=gnuplot.data2
export GPpicwidth=$picwidth
export GPpicheight=$picheight
export GPstart=$start
export GPend=$end
export GPdotsize=$dotsize
export GPtitle="Chromosome $m"
export yminval=$yminv
export ymaxval=$ymaxv
export dotcolour=$dotcol
export dotcolour2=$dotcol2
export GPmarker=$marker


if [ $dotsize == "0.1" ] 
then
    dotstyle=0
else
    dotstyle=7
fi
export GPdotstyle=$dotstyle

if [ $drawgrid == "Y" ]
then
	export GPgrid="1"
else
	export GPgrid="0"
fi


export GPlogscale=$logscale


newinputfile="newinputfile.txt"
newinputfile2="newinputfile2.txt"


### preprocess input file for data files large segments
if [ $stepsize -gt 0 ]
then
	# if begin != end, make entry every x coordinates (so a line is drawn
	
	awk 'BEGIN{
		FS="\t"
		OFS="\t"
	}{
		if(FNR>1){
			i=$2
			j=$3
			while(i<=j){
				print $1,i,i,$4
				i+="'"$stepsize"'"
			}
		}
	}
	END{

	}' $inputfile > $newinputfile

	
	if [ $inputfile2 != "None" ]
	then
		awk 'BEGIN{
			FS="\t"
			OFS="\t"
		}{
			if(FNR>1){
				i=$2
				j=$3
				while(i<=j){
					print $1,i,i,$4
					i+="'"$stepsize"'"
				}
			}
		}
		END{

		}' $inputfile2 > $newinputfile2
		
	fi
else
	cp $inputfile $newinputfile
	if [ $inputfile2 != "None" ]
	then
		cp $inputfile2 $newinputfile2
	else
		newinputfile2="None"
	fi

fi




if [ $chrom == "all" ]
then
	echo "generating plots"
	export GPstart=0
	for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24	
	do
		m=${chrom:3}	
 		end=${chrlen[$m]}
		
		if [ $chrom == "chr23" ]
		then
			chrom="chrX"
			m="X"
		fi
		if [ $chrom == "chr24" ]
		then
			chrom="chrY"
			m="Y"
		fi

		export GPtitle="Chromosome $m"
		export GPend=$end

		###  convert input file
		if [[ $chrprefix == "Y" ]]
		then
			chromval=$chrom
		else
			chromval=${chrom:3}
		fi
		#echo "chromval=$chromval"

		awk 'BEGIN{
				FS="\t"
				OFS="\t"
				count=0	
				binsize="'"$movavg"'"							
			}{				
				if(FNR>1 && $1=="'"$chromval"'")
					print $0
				
			}END{				
			}' $newinputfile > gnuplot.data
			head gnuplot.data

		if [ $inputfile2 != "None" ]
		then
			awk 'BEGIN{
				FS="\t"
				OFS="\t"								
			}{

			
				if(FNR>1 && $1=="'"$chromval"'")
					print $0
			
			}END{				
	
			}' $newinputfile2 > gnuplot.data2
	
			### start GNUplot
			gnuplot ${scriptdir}/scatterplot_2tracks.gnu	
			#gnuplot /data/galaxy-dist/tools/trait/testtool/frequencyplot.gnu	
			wait
			mv out.png chr$m.png
		else
			### start GNUplot
			gnuplot ${scriptdir}/scatterplot.gnu	
			#gnuplot /data/galaxy-dist/tools/trait/testtool/frequencyplot.gnu	
			wait
			mv out.png chr$m.png
		fi

		
	done	

	echo "making montages"
	if [ $subpages == "Y" ]
	then 
		montage chr1.png chr2.png chr3.png chr4.png chr5.png chr6.png chr7.png chr8.png  -geometry ${picwidth}x${picheight}+1+1 -tile 2x4 scatterplot_1-8.png 
		montage chr9.png chr10.png chr11.png chr12.png chr13.png chr14.png chr15.png chr16.png -geometry ${picwidth}x${picheight}+1+1 -tile 2x4 scatterplot_9-16.png
		montage chr17.png chr18.png chr19.png chr20.png chr21.png chr22.png chrX.png chrY.png -geometry ${picwidth}x${picheight}+1+1 -tile 2x4 scatterplot_17-Y.png
	else
		montage chr1.png chr2.png chr3.png chr4.png chr5.png chr6.png chr7.png chr8.png chr9.png chr10.png chr11.png chr12.png chr13.png chr14.png chr15.png chr16.png chr17.png chr18.png chr19.png chr20.png chr21.png chr22.png chrX.png chrY.png -geometry ${picwidth}x${picheight}+1+1 -tile 2x12 scatterplot_all.png
	fi

else  #one single chromosome plot

	echo "extracting data from mastervar file"
	echo "mastervarfile: $mastervarfile"
	echo "chrom: $chrom"
	echo "start: $start"
	echo "end: $end"

	###  convert input file
	if [[ $chrprefix == "Y" ]]
	then
		chromval=$chrom
	else
		chromval=${chrom:3}
	fi
	echo "chromval=$chromval"
	
	awk 'BEGIN{
				FS="\t"
				OFS="\t"
				binsize="'"$movavg"'"+0  	#+0 to force numerical type
				count=0							
			}{
			if(FNR<=5) print $0 > "awk.log"
			if("'"$movavg"'" != 0 ){					
					if($1=="'"$chromval"'"){												
						start_coord[count]=$2
						end_coord[count]=$3
						window[count]=$4
						count++
					}
				}
			else{
				if(FNR>1 && $1=="'"$chrom"'")
					print $0
			}	
			}END{
				# calculate moving average and output to file
				if("'"$movavg"'" != 0){
					#get first average, after that out with the old, in with the new
					print "getting moving averages" > "awk.log"
					print "binsize",binsize > "awk.log"
					print "count",count > "awk.log"
					for (point=0; point<binsize; point++){
						print "k: ",point,"value: "window[point] > "awk.log"
						total+=window[point]
						print "first average",total > "awk.log"						
					}
					for(i=binsize/2; i<(count-(binsize/2)) ;i++){    # for every data point
						total -= window[i - (binsize/2)]
						total += window[i + (binsize/2)]															
						
						print "'"$chromval"'",start_coord[i],end_coord[i],total/binsize	
						print "'"$chromval"'",start_coord[i],end_coord[i],total/binsize	> "awk.log"
					}
					
				}

			}' $newinputfile > gnuplot.data
			head $newinputfile
			head gnuplot.data
			cat awk.log

	if [ $inputfile2 != "None" ]
		then
			awk 'BEGIN{
				FS="\t"
				OFS="\t"	
				binsize="'"$movavg"'"+0  	#+0 to force numerical type
				count=0									
			}{

				if("'"$movavg"'" != 0 ){					
					if($1=="'"$chromval"'"){												
						start_coord[count]=$2
						end_coord[count]=$3
						window[count]=$4
						count++
					}
				}
			else{
				if(FNR>1 && $1=="'"$chrom"'")
					print $0
			}	

			}END{
				# calculate moving average and output to file
				if("'"$movavg"'" != 0){
					#get first average, after that out with the old, in with the new
					print "getting moving averages" > "awk.log"
					print "binsize",binsize > "awk.log"
					print "count",count > "awk.log"
					for (point=0; point<binsize; point++){
						print "k: ",point,"value: "window[point] > "awk.log"
						total+=window[point]
						print "first average",total > "awk.log"						
					}
					for(i=binsize/2; i<(count-(binsize/2)) ;i++){    # for every data point
						total -= window[i - (binsize/2)]
						total += window[i + (binsize/2)]															
						
						print "'"$chromval"'",start_coord[i],end_coord[i],total/binsize	
						print "'"$chromval"'",start_coord[i],end_coord[i],total/binsize	> "awk.log"
					}
					
				}

			}' $newinputfile2 > gnuplot.data2
	
			### start GNUplot
			gnuplot ${scriptdir}/scatterplot_2tracks.gnu	
			#gnuplot /data/galaxy-dist/tools/trait/testtool/frequencyplot.gnu	
			wait
			mv out.png scatterplot.png
		else
			### start GNUplot
			gnuplot ${scriptdir}/scatterplot.gnu	
			#gnuplot /data/galaxy-dist/tools/trait/testtool/frequencyplot.gnu	
			wait
			mv out.png scatterplot.png
		fi
		
	
	
fi




