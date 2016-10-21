#/usr/bin/env bash

# Current Directory
DIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Defaults
config="${DIR}/config.json"
memory=3600
ncores=5

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
        *)
        ;;
    esac
    shift
done

# Compile bsub command
bsub_command="bsub -M${memory} -R 'rusage[mem=${memory}] select[mem>${memory}] span[hosts=1]' -n ${ncores} -o job.log"

echo $bsub_command

# Collect irods
snakemake -s ${DIR}/irods.snake  --configfile=${config} --cluster "bsub -M${memory} -R 'rusage[mem=${memory}] select[mem>${memory}] span[hosts=1]' -n ${ncores} -o job.log" --jobs 16 
