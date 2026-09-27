#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(bambu))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5L)
mode <- args[1]
genome <- args[2]
gtf <- args[3]
cores <- as.integer(args[4])
bams <- args[5:length(args)]
stopifnot(mode %in% c("annotated", "unannotated"), cores >= 1L,
          all(file.exists(c(genome, bams, paste0(bams, ".bai")))))

# Scenario 2 explicitly passes NULL, not an empty or reused annotation set.
annotations <- if (mode == "annotated") prepareAnnotations(gtf) else NULL
set.seed(4262)
dir.create("bambu_output", showWarnings = FALSE)
capture.output(sessionInfo(), file = "bambu_output/sessionInfo.txt")
capture.output(list(mode = mode, genome = genome, annotation = gtf,
                    bams = bams, ncore = cores, NDR = "automatic", seed = 4262),
               file = "bambu_output/parameters.txt")
cat("Bambu mode:", mode, "\nAnnotations supplied:", !is.null(annotations), "\n")
started <- Sys.time()
se <- bambu(reads = bams, annotations = annotations, genome = genome,
            ncore = cores, discovery = TRUE, quant = TRUE)
show(se)
cat("Bambu elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
saveRDS(se, "bambu_output/se.rds")
writeBambuOutput(se, path = "bambu_output")
capture.output(show(se), file = "bambu_output/summarized_experiment.txt")
capture.output(metadata(se)$warnings, file = "bambu_output/bambu_warnings.txt")
# Export an explicit gene-count table in addition to Bambu's standard outputs.
genes <- transcriptToGeneExpression(se)
write.table(as.data.frame(assay(genes, "counts")),
            file = "bambu_output/gene_counts.tsv", sep = "\t", quote = FALSE,
            col.names = NA)
