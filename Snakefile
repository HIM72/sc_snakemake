import os
from utils.parse import read_quants, read_qcs

RESULTS_FOLDER = "/nfs/team205/ge2/felipe"
FASTQ_FOLDER = "data/fastq"
CRAM_FOLDER = "data/cram"
LOG_FOLDER = "logs"
# FASTQ_FOLDER = os.path.join(RESULTS_FOLDER, "test")
LUSTRE = "/lustre/scratch110/sanger/ge2"

with open('items.json') as fh:
    items = json.load(fh)

# LANES = ["1", "2"]
# RUNS = ["20933"]
# SAMPLES = [i for i in range(1, 3)]

SAMPLES = list(set([item["avus"]["tag_index"] for item in items]))
RUNS = list(set([item["avus"]["id_run"] for item in items]))
LANES = list(set([item["avus"]["lane"] for item in items]))

rule all:
    input:
        os.path.join(RESULTS_FOLDER, "results", "tpm.csv"),
        os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon", "multiqc_report.html"),
        os.path.join(RESULTS_FOLDER, "qc", "salmon_qc.csv")

rule multiqc_salmon:
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "merged#{sample}", 
               "quant.genes.sf"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon", "multiqc_report.html")
    log:
        os.path.join(LOG_FOLDER, "multiqc_salmon.log")
    shell:
        "multiqc {salmon_results} -o {out_folder}".format(
            salmon_results=os.path.join(RESULTS_FOLDER, "quant"),
            out_folder = os.path.join(RESULTS_FOLDER, "qc", "multiqc_salmon"))

rule read_salmon_qcs:
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "merged#{sample}", 
               "quant.genes.sf"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_FOLDER, "qc", "salmon_qc.csv"),
    log:
        os.path.join(LOG_FOLDER, "parsing", "salmon_qc.log")
    run:
        qcs = read_qcs(pattern=os.path.join(RESULTS_FOLDER, "quant", "*"))
        qcs.to_csv(output[0])

rule read_reads:
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "merged#{sample}", 
               "quant.genes.sf"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_FOLDER, "results", "reads.csv")
    log:
        os.path.join(LOG_FOLDER, "parsing", "reads.log")
    run:
        results = read_quants(os.path.join(RESULTS_FOLDER, "quant", "*"), 
                              cols=["NumReads"])
        results.to_csv(output[0])

rule read_tpm:
    input:
        expand(os.path.join(RESULTS_FOLDER, "quant", "merged#{sample}", 
               "quant.genes.sf"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_FOLDER, "results", "tpm.csv"),
    log:
        os.path.join(LOG_FOLDER, "parsing", "tpm.log")
    run:
        results = read_quants(pattern=os.path.join(RESULTS_FOLDER, "quant", "*"))
        results.to_csv(output[0])

rule quantify:
    input:
        forward=os.path.join(LUSTRE, "merged", 
                             "merged#{sample}_forward.fastq"),
        reverse=os.path.join(LUSTRE, "merged", 
                             "merged#{sample}_reverse.fastq")
    output:
        os.path.join(RESULTS_FOLDER, "quant", "merged#{sample}", 
                     "quant.genes.sf")
    log:
        lambda wildcards: os.path.join(
            LOG_FOLDER, "salmon", "{sample}.log".format(sample=wildcards.sample))
    shell:
        "salmon quant -i /nfs/team205/.scapi/references/human/salmon_index "
        "-g /nfs/team205/.scapi/references/human/human_gene_map.txt "
        "-l IU -1 {input.forward} -2 {input.reverse} -o "
        "quant/merged#{wildcards.sample}"

rule merge:
    input:
        lambda wildcards: [
            os.path.join(FASTQ_FOLDER, 
                         "{run}_{lane}#{sample}_forward.fastq".format(
                             run=RUNS[0], lane=LANES[0],
                             sample=wildcards.sample)),
            os.path.join(FASTQ_FOLDER, 
                         "{run}_{lane}#{sample}_reverse.fastq".format(
                             run=RUNS[0], lane=LANES[1],
                             sample=wildcards.sample))]
    output:
        temp(os.path.join(LUSTRE, "merged",
                  "merged#{sample}_{direction}.fastq"))
    log:
        lambda wildcards: os.path.join(
            LOG_FOLDER, "merge", "{sample}.log".format(sample=wildcards.sample))
        # os.path.join(LOG_FOLDER, "merge", "{run}_{lane}.log")
    shell:
        "cat {input} > {output}"

rule convert_fastq:
    input:
        lambda wildcards: [
            os.path.join(CRAM_FOLDER, "{run}_{lane}#{sample}.cram".format(
                run=wildcards.run, lane=wildcards.lane, sample=wildcards.sample))]
    output:
        forward=os.path.join(FASTQ_FOLDER, "{run}_{lane}#{sample}_forward.fastq"),
        reverse=os.path.join(FASTQ_FOLDER, "{run}_{lane}#{sample}_reverse.fastq")
    log:
        os.path.join(LOG_FOLDER, "fastq_conversion", "{run}_{lane}#{sample}.log")
    shell:
        "samtools sort -m 10G -n -T %s {input} | "
        "samtools fastq -F 0xB00 -1 {output.forward} -2 {output.reverse} -"

