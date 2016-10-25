# Snakemake Pipelines for Single Cell

## Installation

You must have the following on your path

* Salmon (v0.6.0)
* Samtools
* Python
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

The pipeline bash script can be used as follows (assuming your config is in
~/workflows/config.json):

    pipeline.sh -c ~/workflows/config.json

This will create a workflow with default cluster memory of 3600 and 
5 cores.

Memory can cores can be changed with the `-M` or `-n` flags respectively

