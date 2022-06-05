
# RNA vs. ChIP analysis, [Bragdon *et. al.* 2022](https://doi.org/10.1101/2022.05.22.492993)

[Snakemake](https://snakemake.github.io/) workflow used to compare RNA- and CHIP-seq data for the 2022 publication [*Cooperative assembly confers regulatory specificity and long-term genetic circuit stability*](https://doi.org/10.1101/2022.05.22.492993). For the pipeline with raw data, see the Zenodo archive (coming).

## workflow summary

Matches differential occupancy results from ChIP-seq with differential expression results from RNA-seq by looking for ChIP peaks in a window upstream of transcripts, and does some data visualization comparing ChIP and RNA results from matched peaks/transcripts.
