library(tidyverse)
library(ggrepel)
library(ggpmisc)

main = function(
    input_path = "chip_ZF-cumulative-peaklist-v-rna_verified_genes_plus_Venus_low-affinity-ZF-v-high-affinity-ZF.tsv",
    condition_id = "low-affinity-ZF",
    control_id = "high-affinity-ZF",
    factor = "ZF",
    rna_fdr = 0.1,
    chip_fdr = 0.05,
    output_path = "test.pdf"){

    df = read_tsv(input_path,
                  col_types="ciicicddddddddciicicdddddddddd") %>%
        mutate(rna_significant=(rna_log10_padj > -log10(rna_fdr)),
               chip_significant=(chip_log10_padj > -log10(chip_fdr)))

    plot = ggplot(data = df,
           aes(x=chip_log2FC_enrichment,
               y=rna_log2_foldchange)) +
        geom_hline(yintercept=0, color="grey70", size=0.2) +
        geom_vline(xintercept=0, color="grey70", size=0.2) +
        geom_abline(slope=1, intercept=0, color="grey70", size=0.2) +
        geom_rug(sides="b",
                 aes(color=chip_significant),
                 alpha=0.5,
                 size=0.2) +
        geom_rug(sides="l",
                 aes(color=rna_significant),
                 alpha=0.5,
                 size=0.2) +
        geom_point(aes(color=(rna_significant | chip_significant)),
                   alpha=0.5,
                   shape=16,
                   size=0.8) +
        stat_dens2d_filter(aes(color=(rna_significant | chip_significant),
                               label=ifelse(rna_significant | chip_significant,
                                          rna_name,
                                          NA)),
                           geom="label_repel",
                           keep.number=20,
                           h=0.5,
                           n=500,
                         min.segment.length=0,
                         point.padding=0,
                         box.padding=0,
                         label.padding=0,
                         force=0.1,
                         label.size=NA,
                         fill="#FFFFFF80",
                         size=2,
                         fontface="italic") +
        theme(legend.position="none") +
        scale_x_continuous(name=bquote(.(factor) ~
                                           "ChIP enrichment, log"[2] ~
                                           textstyle(frac(.(condition_id),
                                                .(control_id)))),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(name=bquote(atop("RNA-seq,",
                                            "log"[2] ~
                                                textstyle(frac(.(condition_id),
                                                .(control_id))))),
                           breaks=scales::pretty_breaks(4)) +
        scale_color_manual(values=c("grey60", "#440154FF"),
                           guide=FALSE) +
        theme_light() +
        theme(text = element_text(color="black",
                                  size=8),
              axis.title.y = element_text(angle=0,
                                          vjust=0.5),
              axis.text = element_text(color="black"),
              panel.grid =  element_blank())

    ggsave(output_path,
           plot=plot,
           width=16,
           height=9,
           units="cm")
}

main(input_path = snakemake@input[["table"]],
     condition_id = snakemake@wildcards[["condition"]],
     control_id = snakemake@wildcards[["control"]],
     factor = snakemake@params[["factor"]],
     rna_fdr = snakemake@params[["rna_fdr"]],
     chip_fdr = snakemake@params[["chip_fdr"]],
     output_path = snakemake@output[["svg"]])
