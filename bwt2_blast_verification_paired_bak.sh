#!/bin/bash -ue

sampleId=$1
taxid=$2
pathogen_type=$3
pathogen_db_path=$4
dehost_read1=$5
dehost_read2=$6
min_map_quality=$7
plot_coverage=$8
bin_path=$9


echo "sampleId: ${sampleId}"
echo "taxid: ${taxid}"
echo "pathogen_type: ${pathogen_type}"
echo "pathogen_db_path: ${pathogen_db_path}"
echo "dehost_read1: ${dehost_read1}"
echo "dehost_read2: ${dehost_read2}"
echo "min_map_quality: ${min_map_quality}"
echo "plot_coverage: ${plot_coverage}"
echo "bin_path: ${bin_path}"

# sampleId=SRR12486971
# taxid=10298
# pathogen_type=virus
# pathogen_db_path=/mnt/data/henbio/project/bingyuanti/database/henbio
# dehost_read1=SRR12486971.dehost.R1.fq.gz
# dehost_read2=SRR12486971.dehost.R2.fq.gz
# min_map_quality=30
# plot_coverage=1
# bin_path=/mnt/data/henbio/henbio-wf/PathogenDectectKit/bin


# blast verification
bowtie2 -p 8 --quiet --omit-sec-seq --no-discordant --no-unal \
  -x ${pathogen_db_path}/${pathogen_type}/bowtie2_index/${taxid} \
  -1 ${dehost_read1} \
  -2 ${dehost_read2} | \
  samtools view -@ 8 -f2 -F 256 -q ${min_map_quality} -bS - | samtools sort -@ 8 -o ${sampleId}.${pathogen_type}.${taxid}.bam

# has repeat record
bedtools bamtofastq -i ${sampleId}.${pathogen_type}.${taxid}.bam -fq ${sampleId}.${pathogen_type}.${taxid}.fq
cat ${sampleId}.${pathogen_type}.${taxid}.fq | kz | awk '{if($4 < 0.5) print $1}' | sort | uniq  > ${sampleId}.${pathogen_type}.${taxid}.kz.txt

# remove low complexity
samtools view -h ${sampleId}.${pathogen_type}.${taxid}.bam | grep -v -f ${sampleId}.${pathogen_type}.${taxid}.kz.txt | samtools view -bS - > ${sampleId}.${pathogen_type}.${taxid}.kz.bam
# fix bam file
samtools sort -n ${sampleId}.${pathogen_type}.${taxid}.kz.bam | samtools fixmate -m - ${sampleId}.${pathogen_type}.${taxid}.fixed.bam
# remove dup
samtools sort ${sampleId}.${pathogen_type}.${taxid}.fixed.bam | samtools markdup -r - ${sampleId}.${pathogen_type}.${taxid}.markdup.bam
samtools index ${sampleId}.${pathogen_type}.${taxid}.markdup.bam
samtools sort -n ${sampleId}.${pathogen_type}.${taxid}.markdup.bam | samtools fasta -1 ${sampleId}.${pathogen_type}.${taxid}.R1.fa -2 ${sampleId}.${pathogen_type}.${taxid}.R2.fa -0 /dev/null -s /dev/null -n -

# coverage
read_count=$(($(samtools view ${sampleId}.${pathogen_type}.${taxid}.markdup.bam | wc -l) / 2))
mapped_chrom_count=$(samtools idxstats ${sampleId}.${pathogen_type}.${taxid}.markdup.bam | awk '$3!=0' | wc -l)
samtools depth -a ${sampleId}.${pathogen_type}.${taxid}.markdup.bam > ${sampleId}.${pathogen_type}.${taxid}.depth.txt
seqs=$(samtools idxstats ${sampleId}.${pathogen_type}.${taxid}.markdup.bam | awk '$3!=0' | cut -f1 | head -20 | paste -s -d '|')
gemone_size=$(cat ${pathogen_db_path}/${pathogen_type}/${pathogen_type}.genome.infor.xls | awk -F'\t' -v taxid="${taxid}" '{if($1 == taxid) print $7}')
cat ${sampleId}.${pathogen_type}.${taxid}.depth.txt | awk -v sampleId="${sampleId}" -v taxid="${taxid}" -v seqs="${seqs}" -v mapped_chrom_count="${mapped_chrom_count}" -v read_count="${read_count}" -v gemone_size="${gemone_size}" '{map_bases+=$3; if ($3>0) $3=1; cover_pos+=$3} END \
  {printf "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.4fX\t%.4f\n",sampleId,taxid,seqs,mapped_chrom_count,read_count,gemone_size,map_bases,cover_pos,map_bases/gemone_size,cover_pos/gemone_size*100}' > ${sampleId}.${pathogen_type}.${taxid}.bwt2.txt

# if the bam is null，coverage and dept result will be -nan
sed -i 's/-nan/0/g' ${sampleId}.${pathogen_type}.${taxid}.bwt2.txt

# coverage plot
if [ "${plot_coverage}" -eq 1 ]; then
  for seq in $(samtools idxstats ${sampleId}.${pathogen_type}.${taxid}.markdup.bam | awk '$3!=0' | cut -f1 | head -20)
  do
    samtools depth -a -r ${seq} ${sampleId}.${pathogen_type}.${taxid}.markdup.bam > ${sampleId}.${pathogen_type}.${taxid}_${seq}_depth.txt
  
    Rscript ${bin_path}/coverage_plot.R \
    ${sampleId} \
    ${taxid} \
    ${seq} \
    ${sampleId}.${pathogen_type}.${taxid}_${seq}_depth.txt \
    ${pathogen_type}
  done
fi


# blast verification
blastn -query ${sampleId}.${pathogen_type}.${taxid}.R1.fa \
  -db ${pathogen_db_path}/${pathogen_type}/blast_index/${taxid} \
  -evalue 1e-5 \
  -max_target_seqs 1 \
  -outfmt '6 qseqid sseqid  pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out ${sampleId}.${pathogen_type}.${taxid}.R1.blast.txt
sed -i '1i\qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' ${sampleId}.${pathogen_type}.${taxid}.R1.blast.txt

blastn -query ${sampleId}.${pathogen_type}.${taxid}.R2.fa \
  -db ${pathogen_db_path}/${pathogen_type}/blast_index/${taxid} \
  -evalue 1e-5 \
  -max_target_seqs 1 \
  -outfmt '6 qseqid sseqid  pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out ${sampleId}.${pathogen_type}.${taxid}.R2.blast.txt
sed -i '1i\qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' ${sampleId}.${pathogen_type}.${taxid}.R2.blast.txt


# if R1 or R2 result is null，set it to 0
if [ ! -s "${sampleId}.${pathogen_type}.${taxid}.R1.blast.txt" ] || [ ! -s "${sampleId}.${pathogen_type}.${taxid}.R2.blast.txt" ]; then
  echo 'sample_id  taxid blast_reads avg_pident  min_pident  tmax_pident' | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' > ${sampleId}.${pathogen_type}.${taxid}.blast.res.txt
  echo "${sampleId} ${taxid}  0 0  0 0" | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' >> ${sampleId}.${pathogen_type}.${taxid}.blast.res.txt
else
  Rscript ${bin_path}/blast_result_processed_paired.R \
    ${sampleId} \
    ${taxid} \
    ${sampleId}.${pathogen_type}.${taxid}.R1.blast.txt \
    ${sampleId}.${pathogen_type}.${taxid}.R2.blast.txt \
    ${pathogen_type}
fi
