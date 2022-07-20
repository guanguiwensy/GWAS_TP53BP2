#!/bin/bash

#Mapping
#$ref means bwa reference index; $i means sample ID; $k means lane ID

ref=$3

for i in $(seq $1 $2)
do
  echo $i
        ls S${i}.*1.fq.gz >1
        ls S${i}.*2.fq.gz >2
        paste 1 2 >config
  k=1
        cat config |while read id
        do
                arr=($id)
                fq1=${arr[0]}
                fq2=${arr[1]}
                bwa mem -t 42 -M -R "@RG\tID:S$i_$k\tSM:S$i\tLB:WGS\tPL:Illumina" $ref $fq1 $fq2 >${fq1}.sam
                gatk --java-options "-Xmx60G -Djava.io.tmpdir=./"  SortSam --MAX_RECORDS_IN_RAM 20000000 -SO coordinate -I ${fq1}.sam -O /home/guanguiwen/data2/${fq1}.bam
                rm ${fq1}.sam
    k=$((k+1))
  done
        bam=/home/guanguiwen/data2/S${i}.*.bam
        samtools merge -@44 S${i}.merge.bam $bam
done

#Variant calling

i=*.merge.marked.bam

freebayes-parallel <(fasta_generate_regions.py /home/guanguiwen/data1/reference/human/hg38/Homo_sapiens_assembly38.fasta.fai 100000) 22 \
-f /home/guanguiwen/data1/reference/human/hg38/Homo_sapiens_assembly38.fasta $i >/home/guanguiwen/data2/merge.vcf

