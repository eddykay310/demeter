#!/bin/bash

seqid=$1
seqid_short=${seqid##*/}
organism_tag=$2
id="${organism_tag}${seqid_short:3}"
keep_transposease=$3

# printf "seqid $1\norg tag $2\nkeep $3\n"

mkdir -p "${seqid}/pseudogenome"

# eliminating transposase and duplicates in pseudofinder
ids=$( grep ">" ${seqid}/cdhit/${id}_cdhit.fasta | tr -d '>' | cut -d' ' -f1 | grep "PF" | sed 's/PF_//' )
echo "$(wc -l <<< $ids) PF ids found"
wanted=$( egrep -h "$ids" ${seqid}/rt_output/prokka_annotation/${id}.faa )
grep -v "hypothetical protein\|putative protein" <<< "$wanted" | sort | uniq -cf 1 | sed 's/^ *//g' | sed 's/^ *\([0-9][0-9]*\) *\(.*\)/\2 \(\1\)/' | sed 's/_/_PF_/' > ${seqid}/pseudogenome/pseudogenome_PF_unique.fasta
if [[ "$keep_transposease" = "F" ]]; then 
    grep "hypothetical protein\|putative protein" <<< "$wanted" | sed 's/_/_PF_/' >> ${seqid}/pseudogenome/pseudogenome_PF_unique.fasta
else
    grep "hypothetical protein\|putative protein\|transposease" <<< "$wanted" | sed 's/_/_PF_/' >> ${seqid}/pseudogenome/pseudogenome_PF_unique.fasta
fi
echo "$(grep ">" ${seqid}/pseudogenome/pseudogenome_PF_unique.fasta | wc -l) PF pseuds found from prokka protein fasta"

# eliminating transposase and duplicates in DFAST
ids=$( grep ">" ${seqid}/cdhit/${id}_cdhit.fasta | tr -d '>' | cut -d' ' -f1 | grep "DF")
echo "$(wc -l <<< $ids) DF ids found"
wanted=$( egrep -h "$ids" ${seqid}/rt_output/DFAST_output/protein.faa )
grep -v "hypothetical protein\|putative protein" <<< "$wanted" | sort | uniq -cf 1 | sed 's/^ *//g' | sed 's/^ *\([0-9][0-9]*\) *\(.*\)/\2 \(\1\)/' > ${seqid}/pseudogenome/pseudogenome_DF_unique.fasta
if [[ "$keep_transposease" = "F" ]]; then 
    grep "hypothetical protein\|putative protein" <<< "$wanted" >> ${seqid}/pseudogenome/pseudogenome_DF_unique.fasta
else
    grep "hypothetical protein\|putative protein\|transposease" <<< "$wanted" >> ${seqid}/pseudogenome/pseudogenome_DF_unique.fasta
fi
echo "$(grep ">" ${seqid}/pseudogenome/pseudogenome_DF_unique.fasta | wc -l) DF pseuds found from DFAST protein fasta"

# eliminating transposase and duplicates in prokka
ids=$( grep ">" ${seqid}/cdhit/${id}_cdhit.fasta | tr -d '>' | cut -d' ' -f1 | grep -v "DF" | grep -v "PF" )
echo "$(wc -l <<< $ids) Prokka ids found"
wanted=$( egrep -h "$ids" ${seqid}/rt_output/prokka_annotation/${id}.faa )
grep -v "hypothetical protein\|putative protein" <<< "$wanted" | sort | uniq -cf 1 | sed 's/^ *//g' | sed 's/^ *\([0-9][0-9]*\) *\(.*\)/\2 \(\1\)/' > ${seqid}/pseudogenome/pseudogenome_prok_unique.fasta
if [[ "$keep_transposease" = "F" ]]; then 
    grep "hypothetical protein\|putative protein" <<< "$wanted" >> ${seqid}/pseudogenome/pseudogenome_prok_unique.fasta
else
    grep "hypothetical protein\|putative protein\|transposease" <<< "$wanted" >> ${seqid}/pseudogenome/pseudogenome_prok_unique.fasta
fi
echo "$(grep ">" ${seqid}/pseudogenome/pseudogenome_prok_unique.fasta | wc -l) Prokka pseuds found from prokka protein fasta"

#putting it together
cat ${seqid}/pseudogenome/pseudogenome_PF_unique.fasta ${seqid}/pseudogenome/pseudogenome_prok_unique.fasta ${seqid}/pseudogenome/pseudogenome_DF_unique.fasta > ${seqid}/combined_pseudogenome_unique.txt
#grep -v "hypothetical protein\|putative protein" < ${seqid}/combined_pseudogenome_unique.txt | sort -k2 -o ${seqid}/combined_pseudogenome_unique.txt
#cat ${seqid}/combined_pseudogenome_unique.txt ${seqid}/pseudogenome_DF_unique.fasta > ${seqid}/combined_pseudogenome_unique.txt
echo "Pseudogenome size = $(grep -c ">" ${seqid}/combined_pseudogenome_unique.txt)"

source ./writing_pseudogenome_files.sh $seqid $seqid_short $organism_tag
