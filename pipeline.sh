#/usr/bin/env bash

# Current Directory
# DIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# DIR="$(dirname "$(readlink -f $0)" )"

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

pwd

# Defaults
config="${DIR}/config.json"
cluster_config="${DIR}/cluster.json"
memory=36000
ncores=5
do_preprocessing=true
do_irods=true
skip_task=""

# Keyword arguments
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -c|--config)
            echo "setting config"
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
        --skip)
            skip_task="$2"
            shift
        ;;
        *)
        ;;
    esac
    shift
done

echo "using config ${config}"

# Compile bsub command
bsub_command="bsub -M{cluster.memory} -R 'rusage[mem={cluster.memory}] select[mem>{cluster.memory}] span[hosts=1]' -n {cluster.n} -o job.log"
irods_bsub_cmd="bsub -M16000 -R 'rusage[mem=16000] select[mem>16000] span[hosts=1]' -n 1 -o job.log"

# Collect irods
if [[ "$skip_task" != "irods" ]];
then
echo "Collecting iRODS data"
snakemake \
    -s ${DIR}/irods.snake  \
    --configfile=${config} \
    --latency-wait 100 \
    --cluster "$irods_bsub_cmd" \
    --jobs 16;
fi

# Perform the pre-processing pipeline
if [[ "$skip_task" != "preprocessing" ]];
then
echo "Preprocessing"
snakemake \
    -s ${DIR}/Snakefile \
    --configfile=${config} \
    --latency-wait 100 \
    --cluster "$bsub_command" \
    --cluster-config ${cluster_config} \
    --jobs 100;
fi
