################# content, after detected circRNAs, 
################# use the backsplice site to realign, abundance, fraction, unique, transcripts
################# 1 merge 4 tools
################# 2 build index
################# 3 realign the unmapped reads and get the abundance reads
################# 4 remove duplication reads and get the abundance
################# 5 circRNA backsplice site cis splice reads abundance
################# 6 circRNA backsplice site best linear transcripts and circRNA fasta sequence
################# 7 circRNA transcripts prediction and support paired-end reads
################# 8 linear RNA abundance analysis



#dir of data and software
dir=$DIR

#software
put_together_ws=$DIR/put_together_ws.pl
bedtools=$DIR/bedtools
hisat2=$DIR/
samtools=$DIR/
#data
reference_sequence=$DIR/GRCh38.p5.genome.fa
reference=$DIR/


#add the samples
for sample in T2;do
cd $dir/$sample
mkdir circRNA_validate

################# 1 merge 4 tools
#acfs
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$8"_MEA""\t"$7"\t"$9}' ./$sample""_acfs/circle_candidates_CBR.bed12 > ./circRNA_validate/$sample""_CBR.txt
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$8"_CBR""\t"$7"\t"$9}' ./$sample""_acfs/circle_candidates_MEA.bed12 > ./circRNA_validate/$sample""_MEA.txt
cat ./circRNA_validate/$sample""_MEA.txt ./circRNA_validate/$sample""_CBR.txt > ./circRNA_validate/$sample""_acfs.txt
sed -i '1i\chr\tstart\tend\tsample\tacfs\tstrand' ./circRNA_validate/$sample""_acfs.txt
#find_circ
awk '{if($7>=2) print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' ./$sample""_fc/$sample.circ_candidates.bed > ./circRNA_validate/$sample.find_circ.txt #unique mapped 3
sed -i '1i\chr\tstart\tend\tsample\tfind_circ\tstrand' ./circRNA_validate/$sample.find_circ.txt
#circRNA_finder STAR
awk '{if($5>=2) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $sample.s_filteredJunctions.bed > ./circRNA_validate/$sample.circRNA_finder.txt #unique mapped 3
sed -i '1i\chr\tstart\tend\tsample\tcircRNA_finder\tstrand' ./circRNA_validate/$sample.circRNA_finder.txt
#mapsplice
awk -F '[~\t]' '{ if($1==$2 && $3>$4 && $3<=$4+1000000 && $6>=2 && ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$4-1"\t"$3"\t"$14"\t"$6"\t"$7; if($1==$2 && $3<$4 && $3>=$4-1000000 && $6>=2&& ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$3-1"\t"$4"\t"$14"\t"$6"\t"$7}' ./$sample""_mapsplice/fusions_candidates.txt | sed 's/++/+/g' | sed 's/--/-/g' > ./circRNA_validate/$sample.mapsplice.txt
sed -i '1i\chr\tstart\tend\tsample\tmapsplice\tstrand' ./circRNA_validate/$sample.mapsplice.txt
cat ./circRNA_validate/$sample""_acfs.txt ./circRNA_validate/$sample.find_circ.txt ./circRNA_validate/$sample.mapsplice.txt ./circRNA_validate/$sample.circRNA_finder.txt | awk '!a[$1"\t"$2"\t"$3]++' > ./circRNA_validate/$sample.4_tools.txt
awk '{if(!/start/) print }' ./circRNA_validate/$sample.4_tools.txt > ./circRNA_validate/$sample.4_tools.bed
awk '{i++; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"i }' ./circRNA_validate/$sample.4_tools.bed > ./circRNA_validate/rest.4_tools.bed
cd circRNA_validate
#remove some candidate which are boundary shift
$bedtools intersect -wao -a rest.4_tools.bed -b rest.4_tools.bed | awk '{if($2!=$9 && $3!=$10) print}' | awk '{if(($2-$9)*($2-$9)<100 && ($3-$10)*($3-$10)<100 )print }' > merge1
awk 'ARGIND==1{a[$1"\t"$2"\t"$3]}ARGIND>1&&!($1"\t"$2"\t"$3 in a ){print $0}' merge1 rest.4_tools.bed > rest1.4_tools.bed #find shift 10 bp
awk '{a=$7;if($7>$14)a=$14;b=$14;if($7>$14)b=$7;print $0,a,b}' merge1 | awk '!a[$(NF-1)]++' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > merge2 #remove repeat   
cat merge2 rest1.4_tools.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > kept_4_tools.bed

################# 2 build index
awk '!a[$1"\t"$2"\t"$3]++' kept_4_tools.bed > 4_tools_sort.bed
awk '{if($3-$2<=125) print }' 4_tools_sort.bed > short125.bed
#reference 
$bedtools getfasta -fi $reference_sequence -bed short125.bed -fo circRNA_short.fa -s
#repeat quadrupling
paste circRNA_short.fa circRNA_short.fa circRNA_short.fa circRNA_short.fa > short.fa
#test long
awk '{if($3-$2>125) print }' kept_4_tools.bed > long125.bed
#reference
$bedtools getfasta -fi $reference_sequence -bed long125.bed -fo circRNA_long.fa -s
#repeat twice
paste circRNA_long.fa circRNA_long.fa > long.fa
cat long.fa short.fa > total.fa
mkdir index
mv total.fa ./index/
cd ./index/
sed 's/\t>.*//g' total.fa | sed 's/\t//g' > total.repeat.fa
$hisat2/hisat2-build total.repeat.fa total.repeat
$samtools/samtools faidx total.repeat.fa
cd ..
awk -F '[()\t]' '{left=int($4/2)-1; right=int($4/2)+1; print $1"("$2")""\t"left"\t"right"\t"$4"\t1\t"$2}' ./index/total.repeat.fa.fai > ./index/ref2bp.bed

################# 3 realign the unmapped reads
cd $dir/$sample
awk 'NR%4==1' $sample.unmap.fq > unmapped.title
$bedtools bamtobed -split -bed12 -i $sample.snp.bam > accepted_hits.bed
awk -F '[@_\t/]' 'NR==FNR{a[$2]=$0;next}NR>FNR{if($4 in a)print $0}' unmapped.title accepted_hits.bed > unmapped.mapped.bed
rm accepted_hits.bed
rm unmapped.title
$hisat2/hisat2 -p 8 --mm --dta -x ./circRNA_validate/index/total.repeat -U $sample.unmap.fq -S ./circRNA_validate/$sample.unmap.sam 2> $sample.circmap.rate 
$samtools/samtools sort -@ 8 -o ./circRNA_validate/$sample.unmap.bam ./circRNA_validate/$sample.unmap.sam
$samtools/samtools rmdup -S ./circRNA_validate/$sample.unmap.bam --reference ./circRNA_validate/index/total.repeat ./circRNA_validate/rmdup.bam
$bedtools bamtobed -split -bed12 -i ./circRNA_validate/rmdup.bam > ./circRNA_validate/rmdup.bed
#sed -i 's/__/\//g' ./circRNA_validate/rmdup.bed
$bedtools bamtobed -split -bed12 -i ./circRNA_validate/$sample.unmap.bam > ./circRNA_validate/accepted_hits.bed
sed -i 's/__/\//g' ./circRNA_validate/accepted_hits.bed    
$bedtools intersect -split -wa -wb -abam ./circRNA_validate/$sample.unmap.bam -b ./circRNA_validate/index/ref2bp.bed -bed | awk '!a[$0]++' > circ_splice.reads
sed -i 's/__/\//g' circ_splice.reads   #get whole circRNA splice reads in the splice site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed circ_splice.reads > tmp.pair.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-20 <= $20 && a[3]+20 >= $21 && $4!=$22 ) print }' tmp.pair.circ_splice.reads > tmp.pair_proper.circ_splice.reads
awk 'ARGIND==1{a[$0]}ARGIND>1&&!($0 in a ){print $0}' tmp.pair_proper.circ_splice.reads tmp.pair.circ_splice.reads > tmp.pair_unproper.circ_splice.reads
#count
#both paired reads in circRNA sequence and they are not align to genome
awk -F '[_\t]' 'NR==FNR{a[$1"\t"$4]=$0;next}NR>FNR{if($1"\t"$4 in a)print a[$1"\t"$4]"\t"$0}' circ_splice.reads ./circRNA_validate/accepted_hits.bed | awk -F '[_\t]' '{if($4==$24 && $6!~$26) print }' > tmp.both_circ_splice.reads
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_circ_splice.reads > tmp.both_circ_splice.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.both_circ_splice.reads1 circ_splice.reads > tmp.both_circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_circ_splice.reads2 > $sample.both.txt
sed -i '1i\circrna both' $sample.both.txt
#count total 
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' circ_splice.reads > tmp.circ_validate.reads
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads circ_splice.reads > tmp.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_splice.reads2 > $sample.total.txt
sed -i '1i\circrna '$sample'' $sample.total.txt
#count proper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_proper.circ_splice.reads > tmp.circ_validate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads1 tmp.pair_proper.circ_splice.reads > tmp.circ_splice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_splice.reads21 > $sample.pp.txt
sed -i '1i\circrna pp' $sample.pp.txt
# count unproper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_unproper.circ_splice.reads > tmp.circ_unvalidate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_unvalidate.reads1 tmp.pair_unproper.circ_splice.reads > tmp.circ_unsplice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_unsplice.reads21 > $sample.up.txt
sed -i '1i\circrna up' $sample.up.txt
mkdir tmp
cp $sample.total.txt ./tmp
mv $sample.both.txt ./tmp
mv $sample.up.txt ./tmp
mv $sample.pp.txt ./tmp
perl $put_together_ws ./tmp | sort -r -k2 > $sample.circ.txt
rm tmp.*
rm -r ./tmp

################# 4 remove duplication reads and get the abundance
cd $dir/$sample
$bedtools intersect -split -wa -wb -abam ./circRNA_validate/rmdup.bam -b ./circRNA_validate/index/ref2bp.bed -bed | awk '!a[$0]++' > unique.circ_splice.reads     
sed -i 's/__/\//g' unique.circ_splice.reads   #get whole circRNA splice reads in the splice site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed unique.circ_splice.reads > tmp.pair.unique.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-20 <= $20 && a[3]+20 >= $21 && $4!=$22 ) print }' tmp.pair.unique.circ_splice.reads > tmp.pair_proper.unique.circ_splice.reads
awk 'ARGIND==1{a[$0]}ARGIND>1&&!($0 in a ){print $0}' tmp.pair_proper.unique.circ_splice.reads tmp.pair.unique.circ_splice.reads > tmp.pair_unproper.unique.circ_splice.reads
#count
#both paired reads in circRNA sequence and they are not align to genome
awk -F '[_\t]' 'NR==FNR{a[$1"\t"$4]=$0;next}NR>FNR{if($1"\t"$4 in a)print a[$1"\t"$4]"\t"$0}' unique.circ_splice.reads ./circRNA_validate/rmdup.bed | awk -F '[_\t]' '{if($4==$24 && $6!~$26) print }' > tmp.both_unique.circ_splice.reads
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_unique.circ_splice.reads > tmp.both_unique.circ_splice.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.both_unique.circ_splice.reads1 unique.circ_splice.reads > tmp.both_unique.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_unique.circ_splice.reads2 > $sample.both.txt
sed -i '1i\circrna both' $sample.both.txt
#count total 
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' unique.circ_splice.reads > tmp.circ_validate.reads
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads unique.circ_splice.reads > tmp.unique.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.unique.circ_splice.reads2 > $sample.unique.total.txt
sed -i '1i\circrna '$sample'' $sample.unique.total.txt
#count proper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_proper.unique.circ_splice.reads > tmp.circ_validate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads1 tmp.pair_proper.unique.circ_splice.reads > tmp.unique.circ_splice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.unique.circ_splice.reads21 > $sample.pp.txt
sed -i '1i\circrna pp' $sample.pp.txt
# count unproper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_unproper.unique.circ_splice.reads > tmp.circ_unvalidate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_unvalidate.reads1 tmp.pair_unproper.unique.circ_splice.reads > tmp.circ_unsplice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_unsplice.reads21 > $sample.up.txt
sed -i '1i\circrna up' $sample.up.txt
mkdir tmp
cp $sample.unique.total.txt ./tmp
mv $sample.both.txt ./tmp
mv $sample.up.txt ./tmp
mv $sample.pp.txt ./tmp
perl $put_together_ws ./tmp | sort -r -k2 > $sample.unique.circ.txt
rm tmp.*
rm -r ./tmp

################# 5 circRNA backsplice site cis splice reads abundance
cd $dir/$sample
awk -v sample="$sample" 'BEGIN{print "circrna\t"sample"\tboth\tpp\tup"}NR==1{for(i=0;i++<NF;)a[$i]=i;next}{print $a["circrna"]"\t"$a[sample]"\t"$a["both"]"\t"$a["pp"]"\t"$a["up"]}' $sample.circ.txt > $sample.tmp.circ.txt
cat $sample.tmp.circ.txt | sed 's/(-)/ nega/g' | sed 's/(+)/ posi/g' | awk -F '[ -\t:]' '{print $1,$2,$3,$0}' | sed 's/nega/-/g' | sed 's/posi/+/g' | sed 's/ /\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"("$5")\t"$6"\t"$5}' | sed '1d' > $sample.circ.bed
awk '{print $1,$2-1,$2+1,$4,"up",$6"\n"$1,$3-1,$3+1,$4,"down",$6 }' $sample.circ.bed | sed 's/ /\t/g' > $sample.tmp.2.circ.bed 
$bedtools coverage -counts -abam $sample.snp.rmdup.bam -b $sample.tmp.2.circ.bed > $sample.tmp.2.circ.bedtools.rmdup.depth
awk '{a[$4]=a[$4](a[$4]?"\t":"")$5; b[$4]=b[$4](b[$4]?"\t":"")$7;}END{for (j in a) print j"\t"b[j]}' $sample.tmp.2.circ.bedtools.rmdup.depth | sed '1i\circrna\tleft\tright' > $sample.tmp.3.circ.bedtools.rmdup.depth
awk 'NR==FNR{a[$1]=$2"\t"$3;next}NR>FNR{if($1 in a)print $0"\t"a[$1]}' $sample.tmp.3.circ.bedtools.rmdup.depth $sample.tmp.circ.txt > $sample.circ.fraction

################# 6 circRNA backsplice site best linear transcripts and circRNA  
# get circRNA inner exon sequence
$bedtools closest -s -a $sample.circ.bed -b $reference/gencode.v24.refflat_chr.gtf -d | awk '{ if($NF==0) print }' > tmp.1  # circ and exon same as acfs gtf file, inner exon of circRNAs
awk '{print $7"\t"$10"\t"$11"\t"$15"\t"$12"\t"$13}' tmp.1 | awk '{$2=$2; print }' | sed 's/ /\t/g' > tmp.2 #circ related exon
$bedtools getfasta -name $4 -s -tab -fi $reference_sequence -bed tmp.2 -fo tmp.3 #circ related exon fasta
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($15 in a)print $0"\t"a[$15]}' tmp.3 tmp.1 | sort -n -k17 > tmp.4 #match the exon to related circ acording to the transcript and exon name
sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '{if(($3-$2)>=$9) print }' | awk '!a[$4]++' > $sample.tmp.tab2 #circRNA related exon length  make sure each circRNA are output according this condition
sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '!a[$4]++' > $sample.tmp.tab1   #keep the longest
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.tmp.tab2 $sample.tmp.tab1 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($3-$2)"\t"$10"\t"$11}' > $sample.tmp.tab3   #filter some annotation length longer than circRNA spans
cat $sample.tmp.tab3 $sample.tmp.tab2 > $sample.circ_candidates.tab
sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '{if(($3-$2)>$9) print }' > $sample.circ_candidates_multi.tab
awk 'BEGIN{n=10}{for(i=1;i<n;i++)printf $i"\t";print $i}' $sample.circ_candidates.tab | sed 's/>//g' > tmp.anno #circRNA related annotation
#add the utr information
$bedtools closest -s -a tmp.2 -b $reference/utr.bed -d | awk '{ if($NF==0) print }' | sed 's/NM_/NM/g' > tmp.utr_in_exon # circRNA related exon which located in utr region
sed 's/___/\t/g' tmp.utr_in_exon | awk '{if($4==$17) print }' | awk '!a[$4"\t"$5]++' > tmp.utr_in_exon1
awk 'NR==FNR{a[$4]=$0;next}NR>FNR{if($7 in a)print $0"\t"a[$7]}' tmp.utr_in_exon1 tmp.anno | awk '{if($2 < $20 && $3 < $20 ) next; if($2 > $21 && $3 > $21 ) next; print }' > tmp.utr.anno #utr contain circRNAs
awk '{if($2 > $20 && $3 < $21) print }' tmp.utr.anno > $sample.tmp.circ_in_utr # only in utr region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.tmp.circ_in_utr tmp.utr.anno > $sample.tmp.utr_cds.anno
#class 1, utr, 2, only exon, 3, not in exon (intron, intergenic, antisense)
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' tmp.1 $sample.circ.bed > tmp.not_in_exon.bed 
$bedtools closest -s -a tmp.not_in_exon.bed -b $reference/transcript.bed -d | awk '{ if($NF==0) print }' | awk '!a[$4]++' > $sample.tmp.circ_in_intron
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.tmp.circ_in_intron tmp.not_in_exon.bed > tmp.circ_not_intron 
$bedtools closest -S -a tmp.circ_not_intron -b $reference/transcript.bed -d | awk '{ if($NF==0) print }' > $sample.circ_in_antisense 
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.circ_in_antisense tmp.circ_not_intron > $sample.tmp.circ_in_intergenic
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' tmp.utr.anno tmp.anno > $sample.tmp.only_exon.bed
# genetype only focus mRNA or noncoding RNA, NM or NR
#combine together 
#chr start end name abudance strand unique length splice_length annotation_gene annotation_tranccript type(intergenic antisense intron exon 3utr 5utr)
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\tNA\tNA\tintergenic"}' $sample.tmp.circ_in_intergenic > $sample.tmp.circ_in_intergenic1
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\t"$10"\t"$11"\tantisense"}' $sample.circ_in_antisense | awk '!a[$4]++' > $sample.circ_in_antisense1
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\t"$10"\t"$11"\tintron"}' $sample.tmp.circ_in_intron > $sample.tmp.circ_in_intron1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon intron "; if(m==2) type="exon "; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type""$10}' $sample.tmp.only_exon.bed > $sample.tmp.only_exon.bed1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon intron "; if(m==2) type="exon"; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type" "$26" "$10}' $sample.tmp.circ_in_utr > $sample.tmp.circ_in_utr1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon cds intron"; if(m==2) type="exon cds"; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type" "$26" "$10}' $sample.tmp.utr_cds.anno > $sample.tmp.utr_cds.anno1
cat $sample.tmp.circ_in_intergenic1 $sample.circ_in_antisense1 $sample.tmp.circ_in_intron1 $sample.tmp.only_exon.bed1 $sample.tmp.utr_cds.anno1 $sample.tmp.circ_in_utr1 > $sample.annotation.txt
rm *tmp*

################# 7 circRNA transcripts prediction and support paired-end reads
cd $dir/$sample
#according to the paired end reads, to build the circRNA inner splicing site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed circ_splice.reads > tmp.pair.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-2 <= $20 && a[3]+2 >= $21 && $4!=$22 ) print }' tmp.pair.circ_splice.reads > pair_proper_tmp.circ_splice # get the circRNA paired reads
awk 'BEGIN{n=30}{for(i=1;i<30;i++)printf $i"\t";print $i}' pair_proper_tmp.circ_splice | awk '{$23=$22; $22=$1; print }' | sed 's/ /\t/g' | awk 'BEGIN{n=11}{for(i=(NF-n);i<NF;i++)printf $i"\t";print $i}' > pair_proper_tmp.circ_splice.tmp.1 # convert it to bed12
awk '{if($10==2) {n11=split($11,t11,",");n12=split($12,t12,","); print $1"\t"($7+t11[1])"\t"($8-t11[2])"\t"$4"\t"$5"\t"$6} }' pair_proper_tmp.circ_splice.tmp.1 > pair_proper_tmp.circ_splice.tmp.2 # transform the split information to splice site, 2 split
awk '{if($10==3) {n11=split($11,t11,",");n12=split($12,t12,","); for (i = 1; i < 3; i++) {print $1"\t"($7+t11[i]+t12[i])"\t"($7+t12[i+1])"\t"$4"\t"$5"\t"$6}} }' pair_proper_tmp.circ_splice.tmp.1 > pair_proper_tmp.circ_splice.tmp.3 # transform the split information to splice site, 3 split
#get the circRNA inner splice site
cat pair_proper_tmp.circ_splice.tmp.2 pair_proper_tmp.circ_splice.tmp.3 | awk '!a[$1"\t"$2"\t"$3"\t"$4"\t"$5]++' > pair_proper_tmp.circ_splice.inner  # inner has splice site other there no detected splice site. maybe owning to low abundance
awk '{a[$1"\t"$2"\t"$3"\t"$4]++}END{for(i in a){print i"\t"a[i] | "sort -k 1"}}' pair_proper_tmp.circ_splice.inner | awk '{m=substr($4,length($4)-1,1); print $0"\t"m }' > pair_proper.circ_splice.intro  #inner of circRNA splicing site format same as intron
#test whether the left and right are overlap for each circRNAs positive is origin from right, negative is origin from left
awk '{a[$13]=1; if($24~/\-/ && $26>left[$13] ) left[$13]=$26; if($24~/\+/ && (!length(right[$13]) || $25<right[$13])) right[$13]=$25;}END{for(i in a){print i,right[i],left[i]}}' pair_proper_tmp.circ_splice | awk '{if($2<$3) print }' > pair_proper.circ_splice_spanned
#bulid a transcripts acording to the span circRNAs inner splicing site, detection of splicing site. maybe not all, need high abundance
awk -F '[ \t]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($4 in a)print $0}' pair_proper.circ_splice_spanned pair_proper.circ_splice.intro | sort -n -k3 | sort -n -k2 | awk '{a[$4]=a[$4](a[$4]?",":"")$2; b[$4]=b[$4](b[$4]?",":"")$3; c[$4]=c[$4](c[$4]?",":"")$5; d[$4]+=1; e[$4]+=$5; }END{for (j in a) print j"\t"d[j]"\t"e[j]"\t"c[j]"\t"a[j]"\t"b[j]}' > $sample.circ.transcripts
#compare with intron, extract the intron bed format
#awk '{ n9 = split($9, t9, ",");n10 = split($10, t10, ","); for (i = 0; ++i < n9-1;) { print $2"\t"t10[i]"\t"t9[i + 1]"\t"i "I@" $1"\t"$8"\t"$3 }}' gencode.v24.all.refFlat.txt > gencode.v24.all.refFlat.intron.bed
#awk '{ if($8==1) print $2"\t"$4"\t"$5"\t"$1"\t"$8"\t"$3}' gencode.v24.all.refFlat.txt > gencode.v24.all.refFlat.single.bed
$bedtools intersect -wa -wb -b $reference/gencode.v24.all.refFlat.intron.bed -a pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_intron  #overlap with intron
awk '{if($2==$8 && $3==$9) print }' pair_proper_tmp.circ_splice.in_intron > pair_proper_tmp.circ_splice.in_intron_match
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron_match pair_proper_tmp.circ_splice.in_intron > pair_proper_tmp.circ_splice.in_intron_mis
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_not_intron # not overlap with intron
awk '!a[$4]++' pair_proper_tmp.circ_splice.in_intron_mis > pair_proper.circ_splice.in_intron_mis
#not in intron region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_out
$bedtools intersect -wa -wb -b $reference/gencode.v24.all.refFlat.single.bed -a pair_proper_tmp.circ_splice.in_out > pair_proper.circ_splice.in_single
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' pair_proper.circ_splice.in_single pair_proper_tmp.circ_splice.in_out > pair_proper.circ_splice.intergenic

#extract some discordant alignments
#circRNA bed format
sed 's/(-)/(n)/g' $sample.circ.fraction | sed 's/(+)/(p)/g' | awk '{ print $1,$2"_"$4 }' | awk -F '[ -:\(\)]' '{ print $1,$2,$3,$0,$4}' | sed 's/p/+/g' | sed 's/n/-/g' | sed 's/ /\t/g' | sed '1d' > $sample.tmp.validate.circ.bed
$bedtools intersect -f 0.99 -split -abam $sample.snp.bam -b $sample.tmp.validate.circ.bed > $sample.tmp.validate.circ.bam   #extract the reads in the circRNA region
$samtools/samtools view $sample.tmp.validate.circ.bam | awk '{if(($2==161 && $9<0) || ($2==81 && $9>0) || ($2==97 && $9>0) || ($2==145 && $9<0) ) print }' > $sample.tmp.validate.circ.support #get the support paired circRNAs
awk 'a[$1]++' $sample.tmp.validate.circ.support > $sample.tmp.validate.circ.support.title #filter some one end reads. kept both paired end reads
awk -F '[\t/]' 'NR==FNR{a[$1]=1;next}NR>FNR{if($1 in a)print $0}' $sample.tmp.validate.circ.support.title $sample.tmp.validate.circ.support | sort -k 1 > $sample.tmp.validate.circ.support1 #get the proper paired reads
$samtools/samtools view -H $sample.snp.bam > tmp.head.sam
cat tmp.head.sam $sample.tmp.validate.circ.support1 | $samtools/samtools view -bS > $sample.tmp.validate.circ.support.bam  # return to bam for further analysis bed split
$bedtools bamtobed -split -bed12 -i $sample.tmp.validate.circ.support.bam > $sample.tmp.validate.circ.support.bed
$bedtools intersect -wa -wb -a $sample.tmp.validate.circ.support.bed -b $sample.tmp.validate.circ.bed | awk '{if($2>$14 && $3<$15) print }' > $sample.tmp.validate.circ.match
# match the circRNA with related match paired reads 
awk 'ARGIND==1{a[$1]}ARGIND>1&&!($1 in a ){print $0}' $sample.tmp.validate.circ.match $sample.tmp.validate.circ.support.bed > $sample.tmp.validate.circ.support.test # test it whether it is match with circRNAs it should be no
awk '{if($10==2) {n11=split($11,t11,",");n12=split($12,t12,","); print $1"\t"($7+t11[1])"\t"($8-t11[2])"\t"$16"\t"$4"\t"$6} }' $sample.tmp.validate.circ.match > $sample.tmp.validate.circ.match.2  # transform the reads to splice site
awk '{a[$1"\t"$2"\t"$3"\t"$4]++}END{for(i in a){print i"\t"a[i] | "sort -k 1"}}' $sample.tmp.validate.circ.match.2 | awk '{m=substr($4,length($4)-1,1); print $0"\t"m }' > $sample.validate.circ.match #circRNA related splicing site and support paired reads number
#comapre the support with circRNA other end.
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]=$0}ARGIND>1&&($1"\t"$2"\t"$3"\t"$4 in a ){print $0"\t"a[$1"\t"$2"\t"$3"\t"$4]}' $sample.validate.circ.match pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.paired_support #different
#compare with intron
$bedtools intersect -wa -wb -b $reference/gencode.v24.all.refFlat.intron.bed -a $sample.validate.circ.match > $sample.validate.circ.match_intron
awk '{if($2==$8 && $3==$9) print }' $sample.validate.circ.match_intron > $sample.validate.circ.match_intron_match
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' $sample.validate.circ.match_intron_match $sample.validate.circ.match_intron > $sample.validate.circ.match_intron_mis
awk '!a[$4]++' $sample.validate.circ.match_intron_mis | wc -l
awk '{if($2==$8 || $3==$9) print }' pair_proper_tmp.circ_splice.in_intron_mis | awk '!a[$4]++' | wc -l
#not in intron region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.validate.circ.match_intron $sample.validate.circ.match > $sample.validate.circ.match.in_out
$bedtools intersect -wa -wb -b $reference/gencode.v24.all.refFlat.single.bed -a $sample.validate.circ.match.in_out > $sample.validate.circ.match.in_single
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' $sample.validate.circ.match.in_single $sample.validate.circ.match.in_out > $sample.validate.circ.match.in_intergenic
rm *tmp*

done




