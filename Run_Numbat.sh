conda activate Numbat
pileup_and_phase="pileup_and_phase.R"
gamp="genetic_map_hg38_withX.txt"
snpvcf='genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf'
hg38_1000g="1000G_hg38"
ncores=16
myPath=/media/MaleBRCA
OUTdir=$myPath/result/OUT_Numbat
mkdir -p $OUTdir
cd $OUTdir

if [ -f $myPath/Numbat.list ]; then
  rm $myPath/Numbat.list
  ls -d $myPath/data/*_cellranger710/ >> $myPath/Numbat.list
else
  ls -d $myPath/data/*_cellranger710/ >> $myPath/Numbat.list
fi

cat $myPath/Numbat.list | while read id
do
SAMPLEID=$(echo "$id" | cut -d'/' -f7)
echo $SAMPLEID
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
barcodes=$cellranger_outDir/filtered_feature_bc_matrix/barcodes.tsv.gz
mkdir $SAMPLEID
Rscript $pileup_and_phase \
    --label ${SAMPLEID} \
    --samples ${SAMPLEID} \
    --bams $BAM \
    --barcodes $barcodes \
    --outdir $OUTdir/$SAMPLEID \
    --gmap $gamp \
    --snpvcf $snpvcf \
    --paneldir $hg38_1000g \
    --ncores $ncores 
done