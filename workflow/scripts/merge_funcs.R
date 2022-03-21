#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(stringr)
library(broom)

# Logging file
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

kegg <- read.delim(snakemake@input[["coverage"]])
scaffolds2bin <- read.delim(snakemake@input[["scaffolds"]])
edited_kegg <- inner_join(kegg, scaffolds2bin)
head(edited_kegg)

gtdbtk <- read.delim(snakemake@input[["gtdbtk"]])
head(gtdbtk)

# editing the edited_kegg file
edited_kegg$Sample <- gsub("_contig.*","",edited_kegg$Contig)
# edited_kegg$Bin <- paste(edited_kegg$Sample, edited_kegg$Bin, sep="_")
head(edited_kegg)

# merging the KEGG file wih the gtdbtk file
edited_kegg <- inner_join(edited_kegg, gtdbtk)
write.csv(edited_kegg, snakemake@output[["merged"]], quote = FALSE)
