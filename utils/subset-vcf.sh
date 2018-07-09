#!/bin/bash

### USAGE
# VCF file should be provided in first argument
# Pedigree file should be provided in second argument
# Desired output location should be in final argument
###


VCF=$1
PED=$2
OUTPUT=$3

echo "Pedigree being used... $PED"
echo "Output directed to... $OUTPUT"
echo "VCF being used... $VCF"
echo ""
num_founders_needed=$(cat $PED | awk 'BEGIN { count = 0 } { if ( $3=="0" && $4=="0" ) {  count++ } } END { print count }')
num_founders_provided=$(grep '^#' -v -m 1 $VCF | awk '{ print NF-9 } ')
echo "Number of founders in VCF file: $num_founders_provided"
echo "Number of founders VCF is being subsetted to: $num_founders_needed"

if [[ $num_founders_needed -gt $num_founders_provided ]]; then
    echo "More founders needed than are provided... No need to subset"
    exit
fi

random_founders=$(shuf -i 0-$num_founders_provided -n $num_founders_needed)
for indiv in $random_founders; do
    echo "i$indiv" >> ./temp-subset
done

cat $VCF | tr ' ' '\t' | vcf-subset -c ./temp-subset > $OUTPUT
echo "Completed"
rm ./temp-subset
