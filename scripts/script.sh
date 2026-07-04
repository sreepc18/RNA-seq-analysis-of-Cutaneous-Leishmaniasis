#!/bin/bash

RESULTS=../results
INDEX=../data/Homo_sapiens.GRCh38.cdna.all.index
STUDY_DESIGN=../data/studydesign.txt

mkdir -p $RESULTS/fastqc
mkdir -p $RESULTS/kallisto
mkdir -p $RESULTS/kallisto_log

fastqc ../data/*.fastq.gz -o $RESULTS/fastqc -t 8

for file in ../data/*.fastq.gz
do
    srr=$(basename "$file" | cut -d'_' -f1)

    sample=$(awk -F'\t' -v srr="$srr" '$2 == srr {print $1}' "$STUDY_DESIGN")
    group=$(awk -F'\t' -v srr="$srr" '$2 == srr {print $3}' "$STUDY_DESIGN")

    if [ -z "$sample" ]; then
        echo "WARNING: No match for $srr in study design. Using SRR ID."
        sample=$srr
        group="unknown"
    fi

    LABEL=${group}_${sample}
    OUT=$RESULTS/kallisto/$LABEL
    LOG=$RESULTS/kallisto_log/${LABEL}.log

    kallisto quant -i "$INDEX" -o "$OUT" --single -l 250 -s 30 "$file" -t 8 &> "$LOG"
    echo "Mapping done: $srr → $LABEL"
done

multiqc -d $RESULTS -o $RESULTS/multiqc