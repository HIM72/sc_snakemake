import os
from utils.parse import read_quants, read_qcs
from utils.file_control import create_file_targets

# Define the directories
DATA_FOLDER = config['data_folder']
CRAM_FOLDER = config.get('cram_folder', os.path.join(DATA_FOLDER, 'cram'))
FASTQ_FOLDER = config.get('fastq_folder', os.path.join(DATA_FOLDER, 'fastq'))
RESULTS_FOLDER = config.get('results_folder', 
                            os.path.join(DATA_FOLDER, 'results'))
LOG_FOLDER = config.get('log_folder', os.path.join(DATA_FOLDER, 'logs'))
TEMP_FOLDER = config.get('tmp_folder', os.path.join(DATA_FOLDER, '.tmp'))
FROM = config.get('from', 'cram')
REVERSE = config.get('reverse', 'reverse')
FORWARD = config.get('forward', 'forward')

# This is the temporary directory for faster alignment
LUSTRE = config['lustre_folder']

# Define the target files
merge_mapper = {}
if FROM == 'cram':
    glob_pattern = os.path.join(CRAM_FOLDER, config['pattern'] + '.cram')
elif FROM == 'fastq':
    glob_pattern = os.path.join(FASTQ_FOLDER, config['pattern'] + '.fastq')
glob_variables = glob_wildcards(glob_pattern)

# Create merged targets, if any
original_samples, final_samples, merge_mapper = create_file_targets(
    glob_variables, config)

rule all:
    input:
        os.path.join(RESULTS_FOLDER, "results", "tpm.csv"),
        os.path.join(RESULTS_FOLDER, "results", "counts.csv"),
        os.path.join(RESULTS_FOLDER, "qc", "salmon_qc.csv")
        # Currently we are leaving out fastqc as it has a tendency to crash
        # Something to do with java VM memory
        # os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon", "multiqc_report.html"),
        # os.path.join(RESULTS_FOLDER, "qc", "multiqc_fastqc", "multiqc_report.html"),

rule multiqc_salmon:
    # Combine the quality outputs into a multiqc report
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "{sample}", 
               "quant.genes.sf"), sample=final_samples)
    output:
        os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon", "multiqc_report.html")
    log:
        os.path.join(LOG_FOLDER, "multiqc_salmon.log")
    shell:
        "multiqc {salmon_results} -o {out_folder}".format(
            salmon_results=os.path.join(RESULTS_FOLDER, "quant"),
            out_folder = os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon"))

rule read_salmon_qcs:
    # Combine the salmon QC output into a single table
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "{sample}", 
               "quant.genes.sf"), sample=final_samples)
    output:
        os.path.join(RESULTS_FOLDER, "qc", "salmon_qc.csv"),
    log:
        os.path.join(LOG_FOLDER, "parsing", "salmon_qc.log")
    run:
        qcs = read_qcs(pattern=os.path.join(RESULTS_FOLDER, "quant", "*"))
        qcs.to_csv(output[0])

rule read_reads:
    # Combine the read count into a single table
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "{sample}", 
               "quant.genes.sf"), sample=final_samples)
    output:
        os.path.join(RESULTS_FOLDER, "results", "counts.csv")
    log:
        os.path.join(LOG_FOLDER, "parsing", "reads.log")
    run:
        results = read_quants(os.path.join(RESULTS_FOLDER, "quant", "*"), 
                              col="NumReads")
        results.to_csv(output[0])

rule read_tpm:
    # Combine the computed TPM into a single table
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "{sample}", 
               "quant.genes.sf"), sample=final_samples)
    output:
        os.path.join(RESULTS_FOLDER, "results", "tpm.csv"),
    log:
        os.path.join(LOG_FOLDER, "parsing", "tpm.log")
    run:
        results = read_quants(pattern=os.path.join(RESULTS_FOLDER, "quant", "*"))
        results.to_csv(output[0])

rule quantify:
    input:
        forward=os.path.join(LUSTRE, "{sample}_" + FORWARD + ".fastq"),
        reverse=os.path.join(LUSTRE, "{sample}_" + REVERSE + ".fastq")
    output:
        os.path.join(RESULTS_FOLDER, "quant", "{sample}", "quant.genes.sf")
    log:
        lambda wildcards: os.path.join(
            LOG_FOLDER, "salmon", "{sample}.log".format(sample=wildcards.sample))
    params:
        out_folder=lambda wildcards: os.path.join(RESULTS_FOLDER, "quant",
                                                  wildcards.sample),
        index=config['index'],
        gene_map=config['gene_map']
    shell:
        "if [ $(wc -c < '{input.forward}') -le 1200 ]; then "
        "touch {params.out_folder}/quant.genes.sf; "
        "else "
        "salmon quant -i {params.index} "
        "-g {params.gene_map} "
        "-l IU "
        "-1 {input.forward} -2 {input.reverse} "
        "-o {params.out_folder}; "
        "fi; "

rule copy_unmerged_forward:
    # If we aren't merging anything, copy the files directly to lustre
    input:
        lambda wildcards: os.path.join(
            FASTQ_FOLDER,
            "{sample}_{direction}.fastq".format(sample=wildcards.sample,
                                                direction=FORWARD))
    output:
        temp(os.path.join(LUSTRE, "{sample}_" + FORWARD + ".fastq"))
    shell:
        "cp {input} {output}"

rule copy_unmerged_reverse:
    # If we aren't merging anything, copy the files directly to lustre
    input:
        lambda wildcards: os.path.join(
            FASTQ_FOLDER,
            "{sample}_{direction}.fastq".format(sample=wildcards.sample,
                                                direction=REVERSE))
    output:
        temp(os.path.join(LUSTRE, "{sample}_" + REVERSE + ".fastq"))
    shell:
        "cp {input} {output}"

rule merge_reverse:
    # Merge the reverse-stranded files
    input:
        lambda wildcards: [
            os.path.join(FASTQ_FOLDER, 
                         "{original_sample}_{direction}.fastq".format(
                           original_sample=merge_mapper[wildcards.sample][0],
                           direction=REVERSE)),
            os.path.join(FASTQ_FOLDER, 
                         "{original_sample}_{direction}.fastq".format(
                           original_sample=merge_mapper[wildcards.sample][1],
                           direction=REVERSE))]
    output:
        temp(os.path.join(LUSTRE, "{sample}_" + REVERSE + ".fastq"))
    log:
        lambda wildcards: os.path.join(
            LOG_FOLDER, "merge", "{sample}.log".format(sample=wildcards.sample))
    shell:
        "cat {input} > {output}"

rule merge_forward:
    # Merge the forward-stranded files
    input:
        lambda wildcards: [
            os.path.join(
                FASTQ_FOLDER, "{original_sample}_{direction}.fastq".format(
                    original_sample=merge_mapper[wildcards.sample][0],
                    direction=FORWARD)),
            os.path.join(
                FASTQ_FOLDER, "{original_sample}_{direction}.fastq".format(
                    original_sample=merge_mapper[wildcards.sample][1],
                    direction=FORWARD))]
    output:
        temp(os.path.join(LUSTRE, "{sample}_" + FORWARD + ".fastq"))
    log:
        lambda wildcards: os.path.join(
            LOG_FOLDER, "merge", "{sample}.log".format(sample=wildcards.sample))
    shell:
        "cat {input} > {output}"

rule multiqc_fastqc:
    # Combine the fastqc reports into a single HTML file
    input:
        lambda wildcards: 
            expand(os.path.join(TEMP_FOLDER, "fastqc", 
                                "{original_sample}_{direction}_fastqc.zip"),
                original_sample=original_samples, 
                direction=[FORWARD, REVERSE])
    output:
        expand(os.path.join(RESULTS_FOLDER, "qc", "multiqc_fastqc",
                            "multiqc_report.html"))
    params:
        fastqc_folder=lambda x: os.path.join(TEMP_FOLDER, "fastqc"),
        multiqc_folder=lambda x: os.path.join(RESULTS_FOLDER, "qc",
                                              "multiqc_fastqc")
    shell:
        "multiqc {params.fastqc_folder} -o {params.multiqc_folder}"

rule fastqc:
    # Run fastqc
    input:
        forward=os.path.join(
            FASTQ_FOLDER, "{original_sample}_" + FORWARD + ".fastq"),
        reverse=os.path.join(
            FASTQ_FOLDER, "{original_sample}_" + REVERSE + ".fastq")
    output:
        temp(os.path.join(
            TEMP_FOLDER, "fastqc", 
            "{original_sample}_" + FORWARD + "_fastqc.zip")),
        temp(os.path.join(
            TEMP_FOLDER, "fastqc", 
            "{original_sample}_" + REVERSE + "_fastqc.zip"))
    params:
        out_dir=lambda wildcards: os.path.join(TEMP_FOLDER, "fastqc")
    shell:
        "fastqc -o {params.out_dir} {input.forward} {input.reverse}"

rule convert_fastq:
    # Convert from cram to fastq
    input:
        lambda wildcards: [
            os.path.join(CRAM_FOLDER, "{original_sample}.cram".format(
                original_sample=wildcards.original_sample))]
    output:
        forward=os.path.join(
            FASTQ_FOLDER, "{original_sample}_" + FORWARD + ".fastq"),
        reverse=os.path.join(
            FASTQ_FOLDER, "{original_sample}_" + REVERSE + ".fastq"),
    log:
        os.path.join(LOG_FOLDER, "fastq_conversion", "{original_sample}.log")
    shell:
        "samtools sort -m 10G -n -T %s {input} | "
        "samtools fastq -F 0xB00 -1 {output.forward} -2 {output.reverse} -"

