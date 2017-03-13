
####related software are in the default direction


#########################  annotation DIR is the user software or data custom fold
RAISE=./
acfs_file=./acfs_example.txt
dir=./

#software
hisat2= DIR
circRNA_finder=DIR
samtools=DIR/samtools
find_circ=DIR/
bedtools=DIR/
acfs=DIR/

#reference
bed=DIR/gencode.v24.all.bed
STAR_index=DIR/
hisat2_index=DIR/hg38_snp_tran_oudejans/genome_snp_tran
find_circ_fa=DIR/chr


MAPSPLICE_DIR=DIR/MapSplice-v2.1.8/
REF_GENOME=DIR/mapsplice_index/
BOWTIE_INDEX=DIR/total
GTF=DIR/gencode.v24.all.gtf





for sample in T2;do

cd $dir
mkdir $sample
mv $sample.*.fq ./$sample
cd ./$sample

perl $circRNA_finder/runStar.pl $sample.1.fq $sample.2.fq $STAR_index ./$sample.
perl $circRNA_finder/postProcessStarAlignment.pl ./ ./
rm $sample.Aligned.out.sam

$hisat2/hisat2 -p 8 --mm --dta -x $hisat2_index -1 ./$sample.1.fq -2 ./$sample.2.fq -S $sample.snp.sam 2> $sample.map.rate 

awk '{ if($2==77 || $2==141 || $2==89 || $2==133 || $2==137 || $2==69) print }' $sample.snp.sam > $sample.unmap.sam
awk '{m=2; if($2==77 || $2==69) m=1; if($6~/\*/) print "@"$1"_"$2"_"m"\n"$10"\n\+\n"$11}' $sample.unmap.sam > $sample.unmap.fq

$samtools sort -@ 8 -o $sample.snp.bam $sample.snp.sam
$samtools rmdup $sample.snp.bam $sample.snp.rmdup.bam
head -n 100000 $sample.snp.sam | awk '{if(/^@/) print }' > head.tmp
rm $sample.snp.sam &

awk '{if($6~/\*/) print }' $sample.unmap.sam > $sample.unmapped.sam

cat head.tmp $sample.unmapped.sam > $sample.unmap.head.sam
$samtools view -Sb $sample.unmap.head.sam > $sample.unmap.bam

python $find_circ/unmapped2anchors.py $sample.unmap.bam > unmap.qfa
rm $sample.unmapped.sam &
rm $sample.unmap.head.sam &
rm $sample.unmap.bam &
rm head.tmp

mkdir $sample""_fc
bowtie2 -p 8 --reorder --mm -M 20 --score-min=C,-15,0 -q -x /results/genome_ucsc/gencode/index/GRCh38.p5.genome -U unmap.qfa 2> $sample.bt2_second.log | $find_circ/find_circ.py -G $find_circ_fa -p $sample -s ./$sample""_fc/$sample.sites.log > ./$sample""_fc/$sample.sites.bed 2> ./$sample""_fc/$sample.sites.reads
rm unmap.qfa
cd ./$sample""_fc
grep circ $sample.sites.bed | grep -v chrM | $find_circ/sum.py -2,3 | $find_circ/scorethresh.py -16 1 | $find_circ/scorethresh.py -15 2 | $find_circ/scorethresh.py -14 2 | $find_circ/scorethresh.py 7 2 | $find_circ/scorethresh.py 8,9 35 | $find_circ/scorethresh.py -17 100000 > $sample.circ_candidates.bed
$bedtools/closestBed -d -a $sample.circ_candidates.bed -b $bed | awk '{if($NF==0)print }' | awk '!a[$4]++' | sort -nr -k5 > $sample.circ_candidates.annotation


cd $dir/$sample   
mkdir $sample""_acfs
awk '{if(NR%4==2) print;}' $sample.unmap.fq > ./$sample""_acfs/unmapped.reads
cd ./$sample""_acfs
awk '{a[$1]++}END{for(i in a){print i,a[i]}}' unmapped.reads > unmapped.group.reads #the methods is too simple to stupid
awk '{i++; print ">"i"reads\n"$1}' unmapped.group.reads > unmapped.group.fa
awk '{i++; print ">"i"reads\t"$2}' unmapped.group.reads |sed 's/>//g' > unmapped.group.expr
sed -i '1i\newid\tsample1' unmapped.group.expr
#rm unmapped.fastq
rm unmapped.reads
cp $RAISE/$acfs_file.txt SPEC_acfs_""$sample.txt
sed -i s/sample1/$sample/g SPEC_acfs_""$sample.txt
perl $acfs/ACF_MAKE.pl SPEC_acfs_""$sample.txt SPEC_acfs_""$sample.sh
sh SPEC_acfs_""$sample.sh
rm circle_candidates_CBR.p1.1
rm circle_candidates_CBR.p1.2
rm circle_candidates_MEA.p1.1
rm circle_candidates_MEA.p1.2
rm circle_candidates_CBR.sam
rm circle_candidates_MEA.sam
rm circle_candidates_CBR.CL.sa
rm circle_candidates_CBR.CL.amb
rm circle_candidates_CBR.CL.ann
rm circle_candidates_CBR.CL.bwt
rm circle_candidates_CBR.CL.pac
rm circle_candidates_CBR.CL
rm circle_candidates_CBR.pseudo.exon.fa
rm circle_candidates_CBR.pseudo.gene.fa
rm circle_candidates_MEA.CL.amb
rm circle_candidates_MEA.CL.ann
rm circle_candidates_MEA.CL.bwt
rm circle_candidates_MEA.CL.pac
rm circle_candidates_MEA.CL.sa
rm circle_candidates_CBR.agtf
rm circle_candidates_CBR.agtf_gene
rm circle_candidates_CBR.refFlat
rm Step4_CBR_finished
rm circle_candidates_MEA.CL
rm circle_candidates_MEA.pseudo.exon.fa
rm circle_candidates_MEA.pseudo.gene.fa
rm circle_candidates_MEA.agtf
rm circle_candidates_MEA.agtf_gene
rm circle_candidates_MEA.err
rm circle_candidates_MEA.ext50
rm circle_candidates_MEA.gtf
rm circle_candidates_MEA.gtf.ext50
rm circle_candidates_MEA.refFlat
rm circle_candidates_MEA.gtf_2G
rm Step4_MEA_finished
rm circle_candidates_CBR
rm circle_candidates_discard
rm circle_candidates_MEA
rm circle_candidates
rm Step4_finished
rm unmap.parsed.2pp.S4
rm Step3_finished
rm Step3_MuSeg_finished
rm unmap.parsed.2pp.S3
rm unmap.parsed.segs.S2.matched
rm unmap.parsed.segs.S2.novel
rm unmap.parsed.segs.S2.novel2
rm unmap.parsed.segs.S2.novel.fa
rm unmap.parsed.segs.S1
rm unmap.parsed.segs.S1.bed12
rm unmap.parsed.segs.S1.fa
rm unmap.parsed.segs.S1.ppparsed.loop
rm unmap.parsed.segs.S2
rm unmap.parsed.segs.S2.bed12
rm unmap.parsed.segs.S2.fa
rm Step2_MuSeg_finished
rm unmap.parsed.segs.S1.ppparsed.id
rm unmap.parsed.segs.S1.pm
rm unmap.parsed.segs.S1.pmparsed
rm unmap.parsed.segs.S1.pp
rm unmap.parsed.segs.S1.ppparsed
rm unmap.parsed.segs_withN
rm Step2_finished
rm unmap.parsed.2pp.S2.sum
rm unmap.parsed.2pp.S2
rm unmap.parsed.2pp.S2.0
rm unmap.parsed.2pp.S2.1
rm unmap.parsed.2pp.S2.2
rm unmap.parsed.2pp.S2.3
rm unmap.parsed.2pp.S2.pgap
rm unmap.parsed.2pp.S2_withN
rm unmap.parsed.1
rm unmap.parsed.2pm
rm unmap.parsed.2pp
rm unmap.parsed.2pp.S1
rm unmap.parsed.multi
rm unmap.parsed.segs
rm unmap.parsed.tmp
rm unmap.parsed.UID
rm unmap.parsed.UID.fa
rm unmap.parsed.unmap
rm Step1_finished
rm unmap.parsed.MT
rm unmap.sam
rm SPEC_acfs_""$sample.sh
rm SPEC_acfs_""$sample.txt
rm unmapped.group.expr
rm unmapped.group.fa
rm unmapped.group.reads

cd $dir/$sample

mkdir $sample""_mapsplice

OUTPUT_DIR=./$sample""_mapsplice
READ_FILE_END1=./$sample.unmap.fq

python $MAPSPLICE_DIR/mapsplice.py \
       -1 $READ_FILE_END1 \
       -c $REF_GENOME \
       -x $BOWTIE_INDEX \
	   --qual-scale phred33 \
       -p 8 --bam --fusion-non-canonical --min-fusion-distance 200 \
	   --gene-gtf $GTF \
       -o $OUTPUT_DIR 2>log.txt
	   


done

