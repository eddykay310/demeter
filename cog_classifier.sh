#!/bin/bash

base_dir=$1
seqid=${base_dir##*/}
organism_tag=$2
id="${organism_tag}${seqid:3}"
cpus=$3

# mkdir -p whole_genome/COG/$id
# COGclassifier -i whole_genome/faas/$id.faa -o whole_genome/COG/$id -t 6 \
#     1> "whole_genome/COG/cog_analysis.log" 2>  >> "whole_genome/COG/cog_analysis_error.log"
#     [ -s whole_genome/COG/cog_analysis_error.log ] && cat whole_genome/COG/cog_analysis_error.log || \
#     echo -e "COG analysis completed\n" | tee -a "whole_genome/COG/cog_analysis.log"

mkdir -p pseudogenome/COG/$id
COGclassifier -i pseudogenome/faas/$id.faa -o pseudogenome/COG/$id -t $cpus \
    1> "pseudogenome/COG/cog_analysis.log" 2> "pseudogenome/COG/cog_analysis_error.log"
    [ -s pseudogenome/COG/cog_analysis_error.log ] && cat pseudogenome/COG/cog_analysis_error.log || \
    echo -e "COG analysis completed\n" | tee -a "pseudogenome/COG/cog_analysis.log"