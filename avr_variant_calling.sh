#!/bin/bash
set -euo pipefail
shopt -s nullglob

WORKDIR="/home/2016042547/AVR_Sequencing/2024_AVR" #adjust accordingly
cd "$WORKDIR"

READS_DIR="reads"
REF="ref/ref_combined.fasta"
THREADS=8

GENES=(AvrSr35 AvrSr50)

module load minimap2
module load samtools
module load medaka
module load bcftools

mkdir -p align medaka_var medaka_vcf mpileup_vcf \
         filtered_vcf per_gene tidy_tables \
         tmp_work logs allele_calls coverage_stats

echo "[INFO] Starting robust pipeline"


# Reference indexing

[ ! -f "${REF}.fai" ] && samtools faidx "$REF"
minimap2 -d ref/ref_combined.mmi "$REF"


# STEP 1: Mapping

echo "[STEP 1] Mapping reads"

for fq in "${READS_DIR}"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)

    minimap2 -t "$THREADS" -ax map-ont ref/ref_combined.mmi "$fq" \
        | samtools sort -@ "$THREADS" -o "align/${sample}.bam"

    samtools index "align/${sample}.bam"
done


# STEP 2: Split BAMs by gene

echo "[STEP 2] Splitting BAMs by gene"

for bam in align/*.bam; do

    # skip indexes
    [[ "$bam" == *.bai ]] && continue

    base=$(basename "$bam" .bam)

    # skip already split BAMs
    [[ "$base" == *_* ]] && continue

    echo "[SPLIT] $base"

    for gene in "${GENES[@]}"; do

        outbam="align/${base}_${gene}.bam"

        samtools view -b "$bam" "$gene" > "$outbam"

        samtools index "$outbam"

    done
done


# STEP 3: Coverage stats

echo "[STEP 3] Coverage stats"

OUT="coverage_stats/summary.tsv"
echo -e "Sample\tGene\tReads\tMedianLen\tDepth" > "$OUT"

for bam in align/*.bam; do
    base=$(basename "$bam" .bam)

    if [[ "$base" != *"_"* ]]; then
        continue
    fi

    sample="${base%%_*}"
    gene="${base##*_}"

    reads=$(samtools view -c -F 260 "$bam")

    median=$(samtools view "$bam" \
        | awk '{print length($10)}' \
        | sort -n \
        | awk '
            {a[NR]=$1}
            END{
                if(NR==0){print 0; exit}
                if(NR%2) print a[(NR+1)/2];
                else print (a[NR/2]+a[NR/2+1])/2
            }')

    depth=$(samtools depth "$bam" \
        | awk '{s+=$3;n++} END{if(n) print s/n; else print 0}')

    echo -e "$sample\t$gene\t$reads\t$median\t$depth" >> "$OUT"
done

echo "[STEP 3] Coverage complete -> $OUT"


# STEP 4


echo "[STEP 4] Medaka"

for bam in align/*_*.bam; do
    base=$(basename "$bam" .bam)
    out="medaka_var/$base"
    mkdir -p "$out"

    medaka_variant -i "$bam" -r "$REF" -o "$out" -t "$THREADS"

    if [ -f "$out/medaka.vcf" ]; then
        bgzip -c "$out/medaka.vcf" > "medaka_vcf/${base}.vcf.gz"
        bcftools index -f "medaka_vcf/${base}.vcf.gz"
    fi
done

echo "[STEP 5] mpileup"

for bam in align/*_*.bam; do
    base=$(basename "$bam" .bam)

    bcftools mpileup -Ou -f "$REF" "$bam" \
        -a FORMAT/DP,FORMAT/AD -d 10000 \
    | bcftools call -mv -Oz -o "mpileup_vcf/${base}.vcf.gz"

    bcftools index -f "mpileup_vcf/${base}.vcf.gz"
done


# STEP 6: Filtering mpileup VCFs

echo "[STEP 6] Filtering mpileup VCFs"

for vcf in mpileup_vcf/*.vcf.gz; do

    [ -e "$vcf" ] || continue

    base=$(basename "$vcf" .vcf.gz)

    echo "[FILTER] $base"

    tmp="tmp_work/${base}.filled.vcf.gz"

    variants=$(bcftools view -H "$vcf" | wc -l)

    if [ "$variants" -eq 0 ]; then
        echo "[WARN] Empty VCF: $base"
        continue
    fi

    # Add VAF (variant allele fraction) tag from AD
bcftools +fill-tags "$vcf" \
    -Oz \
    -o "$tmp" \
    -- -t VAF

bcftools index -f "$tmp"

# Filter using VAF instead of AF
bcftools view \
    -i 'QUAL>=5 && FORMAT/VAF>=0.50 && FORMAT/DP>=50' \
    -Oz \
    -o "filtered_vcf/${base}.vcf.gz" \
    "$tmp"

bcftools index -f "filtered_vcf/${base}.vcf.gz"

done


# STEP 7: Genotype tables

echo "[STEP 7] Genotype scaling from allele depth"

for gene in "${GENES[@]}"; do

    raw_out="per_gene/${gene}.tsv"
    tidy_out="tidy_tables/${gene}.tsv"

    echo -e "Sample\tCHROM\tPOS\tREF\tALT\tGT\tDP\tAF" > "$raw_out"

    echo -e "Sample\tCHROM\tPOS\tREF\tALT\tGT\tDP\tAF\tGT_scale" > "$tidy_out"

    for v in filtered_vcf/*_${gene}.vcf.gz; do

        [ -e "$v" ] || continue

        sample=$(basename "$v" .vcf.gz)

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n' "$v" \
        | awk -v S="$sample" 'BEGIN{OFS="\t"}{

            gt=$5
            dp=$6

            split($7,a,",")
            ref=a[1]
            alt=a[2]

            if((ref+alt)>0)
                af=alt/(ref+alt)
            else
                af=0

            # genotype scaling
            if(af>=0.80)
    g=1
else
    g=0.5

            print S,$1,$2,$3,$4,gt,dp,af >> raw
            print S,$1,$2,$3,$4,gt,dp,af,g >> tidy

        }' raw="$raw_out" tidy="$tidy_out"

    done
done

echo "[DONE] Pipeline complete"
