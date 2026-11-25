#!/usr/bin/env Rscript
# Consensus on the v2 fibroblast matrix. Same fix as 64b: multiData() instead of
# an explicit class, which is what made the call in 71 fail.

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 16)

IN  <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_input_v2"
OUT <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_out_v2"
TAG <- "common_only_u5000d1000"

M <- as.matrix(read.table(file.path(IN, sprintf("M_fib_%s.tsv", TAG)), sep = "\t",
                          header = TRUE, row.names = 1, check.names = FALSE))
samples <- read.table(file.path(IN, sprintf("samples_%s.tsv", TAG)), sep = "\t", header = TRUE)
spp <- unique(samples$species)

sets <- list()
for (sp in spp) {
  cols <- grep(paste0("^", sp, "\\|fib\\|"), colnames(M))
  sets[[sp]] <- t(M[, cols, drop = FALSE])
}
multiExpr <- multiData(
  Zebra = sets[["Zebra"]], Quail = sets[["Quail"]], Chicken = sets[["Chicken"]],
  Duck = sets[["Duck"]], Peacock = sets[["Peacock"]], Turkey = sets[["Turkey"]],
  Emu = sets[["Emu"]]
)
cs <- checkSets(multiExpr)
cat(sprintf("наборов %d, генов %d, образцов по наборам: %s\n",
            cs$nSets, cs$nGenes, paste(cs$nSamples, collapse = ", ")))

gsg <- goodSamplesGenesMS(multiExpr, verbose = 0)
if (!gsg$allOK) {
  cat(sprintf("отброшено генов %d\n", sum(!gsg$goodGenes)))
  for (i in seq_along(multiExpr)) {
    multiExpr[[i]]$data <- multiExpr[[i]]$data[gsg$goodSamples[[i]], gsg$goodGenes, drop = FALSE]
  }
}

pw <- 12
cat(sprintf("беру softPower = %d, как в сети fib скрипта 71\n", pw))

cons <- blockwiseConsensusModules(
  multiExpr, power = pw, networkType = "signed hybrid", TOMType = "signed",
  minModuleSize = 20, deepSplit = 2, mergeCutHeight = 0.25,
  numericLabels = TRUE, maxBlockSize = 20000, checkMissingData = FALSE, verbose = 0
)
cc <- labels2colors(cons$colors)
genes <- colnames(multiExpr[[1]]$data)
write.table(data.frame(gene = genes, module = cc, label_num = cons$colors),
            file.path(OUT, "consensus_modules_fib.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
sz <- as.data.frame(table(cc), stringsAsFactors = FALSE)
names(sz) <- c("module", "n_genes")
write.table(sz, file.path(OUT, "consensus_module_sizes.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat(sprintf("consensus: модулей %d (без grey %d), grey %d генов\n",
            nrow(sz), sum(sz$module != "grey"),
            ifelse("grey" %in% sz$module, sz$n_genes[sz$module == "grey"], 0)))
cat("готово\n")
