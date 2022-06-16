#!/bin/bash
start=`date +%s`

conda activate "bactopia_manually"

help()
{
    echo "Usage: source main.sh [ -f | --fasta ID ]
               [ -rt | --runtype ]
               [ -c | --cpus ]
               [ -rm | --ram ]
               [ -aq | --assembly_qc ]
               [ -ref | --reference_genome ]
               [ -anc | --ancestor_genome ]
               [ -d | --dnds ]
               [ -ot | --organism_tag ]
               [ -kt | --keep_transposease ]
               [ -ps | --pipeline_stage ]
               [ -s | --scaffold_genome ]
               [ -pt | --pseudogene_predictor ]
               [ -h | --help ]"

    set -e
}

SHORT=-f:,rt:,c:,rm:,aq:,ref:,anc:,d:,ot:,kt:,ps:,s:,pt:,h
LONG=fastaID:,runtype:,cpus:,ram:,assembly_qc:,reference_genome:,ancestor_genome:,\
    dnds:,organism_tag:,keep_transposease:,pipeline_stage:,scaffold_genome:,pseudogene_predictor:,help
OPTS=$(getopt -a -n main --options $SHORT --longoptions $LONG -- "$@")

VALID_ARGUMENTS=$# # Returns the count of arguments that are in short or long options

if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  help
fi

eval set -- "$OPTS"

while :
do
  case "$1" in
    -f | --fastaID )
      seqid="$2"
      shift 2
      ;;
    -rt | --runtype )
      runtype="$2"
      shift 2
      ;;
    -c | --cpus )
        cpus="$2"
        shift 2
        ;;
    -rm | --ram )
        ram="$2"
        shift 2
        ;;
    -aq | --assembly_qc )
        assembly_qc_method="$2"
        shift 2
        ;;
    -ref | --reference_genome )
        reference_genome="$2"
        shift 2
        ;;
    -anc | --ancestor_genome )
        ancestor_genome="$2"
        shift 2
        ;;
    -d | --dnds )
        dnds="$2"
        shift 2
        ;;    
    -ot | --organism_tag  )
        organism_tag="$2"
        shift 2
        ;;    
    -kt | --keep_transposease )
        keep_transposease="$2"
        shift 2
        ;;    
    -ps | --pipeline_stage )
        stage="$2"
        shift 2
        ;;    
    -s | --scaffold_genome )
        scaffold="$2"
        shift 2
        ;;   
    -pt | --pseudogene_predictor )
        predictor_type="$2"
        shift 2
        ;;   
    -h | --help)
      help
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      help
      ;;
  esac
done

app_dir=$(realpath main.sh)
mkdir -p $seqid
base_dir=$(realpath $seqid)
echo -e "Starting $seqid\n"
# base_dir=$1
# runtype=$2 #either blank for PE or any letter for SE; preferably SE
[ "$runtype" = "SE" ] && echo "Using single-end reads" || echo "Using paired-end reads"
# cpus=$3
# ram=$4
# assembly_qc_method=$5
! [[ -z "$reference_genome" ]] && reference_genome=$(realpath $reference_genome) || reference_genome="''"
! [[ -z "$ancestor_genome" ]]  &&  ancestor_genome=$(realpath $ancestor_genome) || ancestor_genome="''"
# dnds=$8
# organism_tag=$9
# keep_transposease="${10}"
# stage="${11}"
# scaffold="${12}"

printf "app_dir=app_dir
base_dir=$base_dir
runtype=$2
cpus=$3
ram=$4
assembly_qc_method=$5
reference_genome=$reference_genome
ancestor_genome=$ancestor_genome
dnds=$8
organism_tag=$9
keep_transposease="${10}"
stage="${11}"
scaffold="${12}"
predictor_type="${predictor_type}"\n\n"

re='^[0-9]+$'
if [ $cpus ] ;then cpus=$cpus; elif [[ "$cpus" =~ "$re" ]] ; then echo -e 'Enter a number as input for the number of cpus to be used'; fi

if [ $stage -le 1 ]; then
    echo -e "###########################################################################\n"
    #echo -e "${base_dir##*/} started at $(( (start / 3600) % 60 )):$(( (start % 3600) / 60 )):$(( (start % 3600) % 60 ))" 
    echo -e "${base_dir##*/} started at $(date | awk '{print $4}')"
    echo -e "Gathering reads\n" 
    start_gather=`date +%s`
    ./gather_samples.sh $base_dir
    end_gather=`date +%s`
    runtime=$((end_gather-start_gather))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Gathering reads runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo -e ">$base_dir\nGathering reads,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 2 ]; then
    echo -e "###########################################################################\n"
    echo -e "Reads QC\n"
    start_qc=`date +%s`
    ./qc_reads.sh $base_dir $runtype $cpus
    end_qc=`date +%s`
    runtime=$((end_qc-start_qc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Reads QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Reads QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 3 ]; then
    echo -e "###########################################################################\n"
    echo -e "Assembling genome\n"
    start_assembly=`date +%s`
    ./genome_assembler.sh $base_dir $runtype $cpus $ram $organism_tag
    end_assembly=`date +%s`
    runtime=$((end_assembly-start_assembly))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Genome assembly runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Genome assembly,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 4 ]; then
    echo -e "###########################################################################\n"
    echo -e "Genome assembly QC\n"
    start_assembly_qc=`date +%s`
    # ./assembly_qc_scale.sh $base_dir $cpus $assembly_qc_method $organism_tag
    ./assembly_qc.sh $base_dir $cpus $assembly_qc_method $organism_tag
    end_assembly_qc=`date +%s`
    runtime=$((end_assembly_qc-start_assembly_qc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Genome assembly QC runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Genome assembly QC,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 5 ]; then
    echo -e "###########################################################################\n"
    echo -e "Predicting pseudogenome\n"
    start_pseudogene_proc=`date +%s`
    source ./pseudo_proccessor_main.sh $base_dir $cpus $reference_genome $ancestor_genome $organism_tag $dnds $scaffold $predictor_type
    end_pseudogene_proc=`date +%s`
    runtime=$((end_pseudogene_proc-start_pseudogene_proc))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Pseudogenome processing runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Pseudogenome processing,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 6 ]; then
    conda activate 'bactopia_manually'
    echo -e "###########################################################################\n"
    echo -e "Cleaning and writing pseudogenome\n"
    start_pseudogene_cleaning=`date +%s`
    source ./cleaning_pseudogenome.sh $base_dir $organism_tag $keep_transposease $predictor_type
    end_pseudogene_cleaning=`date +%s`
    runtime=$((end_pseudogene_cleaning-start_pseudogene_cleaning))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "Pseudogenome cleaning runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "Pseudogenome cleaning,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

if [ $stage -le 7 ]; then
    conda activate 'cog_classifier'
    echo -e "###########################################################################\n"
    echo -e "Running COG analysis\n"
    start_COG_analysis=`date +%s`
    source ./cog_classifier.sh $base_dir $organism_tag $cpus
    end_COG_analysis=`date +%s`
    runtime=$((end_COG_analysis-start_COG_analysis))
    hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
    echo "COG analysis runtime: $hours:$minutes:$seconds (hh:mm:ss)"
    echo "COG analysis,$hours:$minutes:$seconds" >> $base_dir/timer.txt
fi

# find $base_dir -maxdepth 1 -type f -name "*.txt" -o -name "*.log" -o -name "total_contigs*" -o -name "*gz" | xargs -I{} rm "{}" 
# find $base_dir -maxdepth 1 -name "shovill" | xargs -I{} rm -r "{}" 
# rm -rf ./test.msh ./*.bt2 ./total_contigs*

# mv -t $base_dir ./test.msh ./*.bt2 ./total_contigs*
# rm -rf $base_dir/*.bt2 $base_dir/*.fastq.gz
mkdir -p "$base_dir/log_files"
mv -t $base_dir/log_files $base_dir/*out.txt $base_dir/*err.txt $base_dir/*processing.txt $base_dir/timer.txt 

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
echo "Pipeline runtime: $hours:$minutes:$seconds (hh:mm:ss)"
echo "$base_dir,$hours:$minutes:$seconds" | tee -a $base_dir/timer.txt timer.txt
echo -e "################################################################################################\n"