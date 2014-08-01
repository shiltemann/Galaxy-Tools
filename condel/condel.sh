#!/bin/bash


#condel.sh $infile $col_sift $col_polyphen $col_mutass $outfile $condel2script $condel3script $condelconfigdir
echo "input parameters: $@"
if [ $# -ne 9 ]
then
	echo "unexpected number of arguments in $0, exiting"
	exit
fi

infile=$1
sift_col=$2
pph2_col=$3
mutass_col=$4
outfile=$5
runMutass=$4

condel2=$6
condel3=$7
condelconfigdir=$8
cgatools=$9



#create version of infile with linenumbers as variant IDs (so we can join file back together later)
awk 'BEGIN{
		FS="\t"
		OFS="\t"
	}{
		if(FNR==1)
			print "tmpvarId", $0
		else
			print FNR,$0

}END{}' $infile > infile_numbered

#added a column, so update sift, pph2 and mutass column variables
id_col=1
sift_col=$[$sift_col+1]
pph2_col=$[$pph2_col+1]
mutass_col=$[$mutass_col+1]


if [[ $runMutass != "None" ]]
then
	echo "running 3-score condel"
	### 3-score Condel ###
	# create input file for condel:  <variantID> <SIFT score> <Polyphen2 score>
	awk 'BEGIN{
	 		FS="\t";
			OFS="\t"
		}{
			if(FNR>1 && $"'"$sift_col"'" != "" && $"'"$pph2_col"'" != "" && $"'"$mutass_col"'" != ""){
				print $"'"$id_col"'",$"'"$sift_col"'",$"'"$pph2_col"'", $"'"$mutass_col"'"
			}	

		}END{}' infile_numbered > condel3in

		# run condel
		echo "tmpvarId	CONDEL_score_SPM	CONDEL_pred_SPM" > condelout
		perl $condel3 $condelconfigdir condel3in >> condelout

		
		#join condel scores back with original file
		echo "joining condel scores with original files"
		$cgatools join --beta \
			--input infile_numbered condelout \
			--output finalout \
			--match tmpvarId:tmpvarId \
			--select A.*,B.CONDEL_score_SPM,B.CONDEL_pred_SPM \
			-a \
			-m compact

else
	echo "running 2-score condel"
	
	### 2-score Condel ###
	# create input file for condel:  <variantID> <SIFT score> <Polyphen2 score>
	awk 'BEGIN{
	 		FS="\t";
			OFS="\t"
		}{
			if(FNR>1 && $"'"$sift_col"'" != "" && $"'"$pph2_col"'" != "" ){
				print $"'"$id_col"'",$"'"$sift_col"'",$"'"$pph2_col"'"
			}	
		}END{}' infile_numbered > condel2in

	# run condel
	echo "tmpvarId	CONDEL_score_SP	CONDEL_pred_SP" > condelout
	perl $condel2 $condelconfigdir condel2in >> condelout
	
	#join condel scores back with original file
	echo "joining condel scores with original file"
	$cgatools join --beta \
		--input infile_numbered condelout \
		--output finalout \
		--match tmpvarId:tmpvarId \
		--select A.*,B.CONDEL_score_SP,B.CONDEL_pred_SP \
		-a \
		-m compact	

fi


#remove first column of output file (our variantID)
cut -f2- finalout > $outfile




