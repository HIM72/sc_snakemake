#/usr/bin/env bash

# Current Directory
DIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Defaults
config="${DIR}/config.json"
memory=36000
ncores=5
do_preprocessing=true
do_irods=true

# Keyword arguments
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -c|--config)
            config="$2"
            shift
        ;;
        -M|--memory)
            memory="$2"
            shift
        ;;
        -n|--ncores)
            ncores="$2"
            shift
        ;;
        -sp|--skip-preprocessing)
            do_preprocessing=false
            shift
        ;;
        -si|--skip-irods)
            do_irods=false
            shift
        ;;
        *)
        ;;
    esac
    shift
done

# Compile bsub command
bsub_command="bsub -M${memory} -R 'rusage[mem=${memory}] select[mem>${memory}] span[hosts=1]' -n ${ncores} -o job.log"

# Collect irods
if $do_irods;
then
snakemake \
    -s ${DIR}/irods.snake  \
    --configfile=${config} \
    --latency-wait 15 \
    --cluster "$bsub_command" \
    --jobs 16;
fi

# Perform the pre-processing pipeline
if $do_preprocessing;
then
snakemake \
    -s ${DIR}/Snakefile \
    --configfile=${config} \
    --latency-wait 15 \
    --cluster "$bsub_command" \
    --jobs 100;
fi
