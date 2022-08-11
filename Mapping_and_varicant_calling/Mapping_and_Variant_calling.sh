#!/bin/bash

#Mapping
func(){
    echo "Usage:"
	echo "Mapping_and_variant_calling.sh [-t threads] [-i index of bwa]" 
	echo "[-v reference vcf] [-f reference genome] [-a reference genome's index with fai format]"
    exit 1
}

while getopts 'h:t:i:v:f:a:' OPT;do
    case $OPT in
	t) threads="$OPTARG";;
	i) index="$OPTARG";;
	v) vcf="$OPTARG";;
	f) fasta="$OPTARG";;
	a) fai="$OPTARG";;
	h) func;;
	?) func;;
	esac
done

for i in $(seq 1 105)
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
                bwa mem -t $threads -M -R "@RG\tID:S$i_$k\tSM:S$i\tLB:WGS\tPL:Illumina" $index $fq1 $fq2 >${fq1}.sam
                gatk --java-options "-Xmx60G -Djava.io.tmpdir=./"  SortSam \
				--MAX_RECORDS_IN_RAM 20000000 -SO coordinate -I ${fq1}.sam -O ${fq1}.bam
                rm ${fq1}.sam
    k=$((k+1))
  done
        bam=S${i}.*.bam
        samtools merge -@ $threads S${i}.merge.bam $bam
		gatk --java-options "-Xmx60G -Djava.io.tmpdir=./" MarkDuplicates \
		--MAX_RECORDS_IN_RAM 20000000 -I S${i}.merge.bam -O S${i}.merge.bam.marked.bam \
        -M S${i}.merge.bam.marked.bam.metrics
    samtools index -@ $threads S${i}.merge.bam.marked.bam
done

#Variant calling

i=*.merge.marked.bam

freebayes-parallel <(fasta_generate_regions.py $fai --chunks 100) $threads \
-g 20000 -f $fasta --min-base-quality 20 --min-mapping-quality 10 -@ $vcf \
-l --use-best-n-alleles 3 $i\
>merge.vcf

bgzip -@ $threads merge.vcf && bcftools index -t --threads $threads merge.vcf.gz;done

