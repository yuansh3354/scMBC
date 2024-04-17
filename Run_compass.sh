conda activate Compass

cat cells | while read filename
do
DIR="${filename%.*}"
meta=meta-${DIR}.csv

compass --data $filename \
	--num-processes 32 \
	--species homo_sapiens \
	--calc-metabolites \
	--num-threads 32 \
	--output-dir $DIR \
	--temp-dir $DIR \
	--lambda 0.25 \
	--microcluster-size 50
done




