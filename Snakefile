#!/usr/bin/env python
import itertools

configfile: "config.yaml"

subworkflow rna_pipe:
    workdir: config["rnaseq_workflow"]

subworkflow chip_pipe:
    workdir: config["chip_workflow"]

CONTROLS = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
CONDITIONS = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

configfile: rna_pipe("config.yaml")
configfile: chip_pipe("config.yaml")

CHIP_ANNOTATIONS = ["peaks"] + list(config["differential_occupancy"]["annotations"].keys() if config["differential_occupancy"]["annotations"] else [])
RNA_ANNOTATIONS = ["transcripts"] + list(config["differential_expression"]["annotations"].keys() if config["differential_expression"]["annotations"] else [])

FACTOR = config["factor"]

localrules:
    target,

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule target:
    input:
        "config.yaml",
        expand(expand("chip_matched_to_rna/chip_{{chip_annotation}}-v-rna_{{rna_annotation}}/{condition}-v-{control}/chip_{{chip_annotation}}-v-rna_{{rna_annotation}}_{condition}-v-{control}.tsv", zip, condition=CONDITIONS, control=CONTROLS), chip_annotation=CHIP_ANNOTATIONS, rna_annotation=RNA_ANNOTATIONS),
        expand(expand("chip_matched_to_rna/chip_{{chip_annotation}}-v-rna_{{rna_annotation}}/{condition}-v-{control}/chip_{{chip_annotation}}-v-rna_{{rna_annotation}}_{condition}-v-{control}.svg", zip, condition=CONDITIONS, control=CONTROLS), chip_annotation=CHIP_ANNOTATIONS, rna_annotation=RNA_ANNOTATIONS)

rule match_chip_to_rna:
    input:
        rna_results = rna_pipe("diff_exp/{rna_annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{rna_annotation}-diffexp-results-all.tsv"),
        chip_results = chip_pipe(f"diff_binding/{{chip_annotation}}/{{condition}}-v-{{control}}/libsizenorm/{{condition}}-v-{{control}}_{FACTOR}-chipseq-libsizenorm-{{chip_annotation}}-diffbind-results-all.tsv"),
        fasta = config["genome"]["fasta"]
    output:
        temp = temp("chip_matched_to_rna/chip_{chip_annotation}-v-rna_{rna_annotation}/{condition}-v-{control}/.chip_{chip_annotation}-v-rna_{rna_annotation}_{condition}-v-{control}.temp"),
        table = "chip_matched_to_rna/chip_{chip_annotation}-v-rna_{rna_annotation}/{condition}-v-{control}/chip_{chip_annotation}-v-rna_{rna_annotation}_{condition}-v-{control}.tsv"
    params:
        upstream = config["search_distance"]["upstream"],
        downstream = config["search_distance"]["downstream"]
    shell: """
        cut -f1-6 {input.rna_results} | \
        tail -n +2 | \
        awk 'BEGIN{{FS=OFS="\t"}} {{($6 == "+") ? $3=$2 + 1 : $2=$3-1}} {{print}}' | \
        bedtools slop -l {params.upstream} -r {params.downstream} -s -i stdin -g <(faidx -i chromsizes {input.fasta}) | \
        paste - <(tail -n +2 {input.rna_results}) > {output.temp}

        bedtools intersect -loj -a {output.temp} -b <(tail -n +2 {input.chip_results}) | \
        cut -f1-6 --complement | \
        cat - <(bedtools intersect -v -a <(tail -n +2 {input.chip_results}) -b {output.temp} | \
            awk 'BEGIN{{FS=OFS="\t"}} {{print ".",-1,-1,".",".",".",".",".",".",".",".",".",".",".",$0}}') | \
        cat <(paste <(head -n 1 {input.rna_results} | sed 's/\</rna_/g') <(head -n 1 {input.chip_results} | sed 's/\</chip_/g')) - > {output.table}
        """

rule datavis:
    input:
        table = "chip_matched_to_rna/chip_{chip_annotation}-v-rna_{rna_annotation}/{condition}-v-{control}/chip_{chip_annotation}-v-rna_{rna_annotation}_{condition}-v-{control}.tsv"
    output:
        svg = "chip_matched_to_rna/chip_{chip_annotation}-v-rna_{rna_annotation}/{condition}-v-{control}/chip_{chip_annotation}-v-rna_{rna_annotation}_{condition}-v-{control}.svg"
    params:
        factor = config["factor"],
        rna_fdr = config["differential_expression"]["fdr"],
        chip_fdr = config["differential_occupancy"]["fdr"]
    conda:
        "envs/plot.yaml"
    script:
        "scripts/rna_v_chip.R"

