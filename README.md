# Snakemake Pipelines for Single Cell

## Installation

This pipeline is based on [snakemake](https://snakemake.readthedocs.io). It is
probably worth reading through the introductory material before using this code.

To get started, you must have the following on your path

* Salmon (ensure that the index specified in your config is the same as the
  salmon in your path)
* Samtools
* Python (snakemake is only supported by Python 3)
* [baton](https://github.com/wtsi-npg/baton)

Make sure the python requirements are fulfilled using the following:

    pip install -r requirements.txt

## How to Use

Generally, full use of the pipeline involves downloading the data from iRODS 
first before processing the data.

These two steps are not done withing the same snakemake workflow,  as the iRODS 
download is limited to 16 simultaneous processes, and in snakemake you cannot 
specify a job limit per step within the same workflow.

Both steps are directed by a shared config.json file. An example file is
included in the repository. Please copy and create your own file per new 
processing task instead of using the provided example file, as this will 
make updating this repository easier.

There are two ways of running the pipeline - either using the bash scripts,
which call snakemake using the required cluster syntax, or by using the
individual .snake files as specified in the original
[repository](https://bitbucket.org/snakemake/snakemake).

### Running the pipeline

The pipeline bash script can be used as follows (assuming your config is in
~/workflows/config.json):

    pipeline.sh -c ~/workflows/config.json


#### The config file

The config file should specify, as a minimum:
* The directory where the cram / fastq files are for analysis
* The lustre folder for temporary files
* The salmon index directory (created with the same version as the salmon in your path)
* A file mapping transcript to gene ID's
* The study ID
* The query of what runs / lanes to fetch from the study
* The file pattern (e.g. `{run}_{lane}#{tag_index}`)

The `query` parameter specifies the iRODS query. The syntax is in JSON.

For example, to get cell with tag index 2 from lane 1 in run 277737, in study 
4132 the `query` parameter should be:

    query: {
        277737: {
            1: [2]
        }
    }

To get all data from lane 1:

    query: {
        277737: {
            1: []
        }
    }

To get all lanes from run 277737:

    query: {
        277737: {}
    }

To get all data from the study, leave the query parameter as an empty object: 

    query: {}

#### Running only part of the pipeline

The two parts of the pipeline can be run independently using the `--skip`
option. To run the iRODS pipeline only, use the following:

    pipeline.sh --skip preprocessing

To skip the iRODS section and preprocess existing files, use the following:

    pipeline.sh --skip irods

The cluster config in this repository adapts to give salmon processes increased
memory and cores.

