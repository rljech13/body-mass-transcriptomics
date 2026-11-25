#!/usr/bin/env Rscript
# Consensus modules over the seven species as seven sets, inside fibroblasts.
#
# The call in script 64 failed with "data does not appear to have the correct
# format" because the multiExpr list was given an explicit class. WGCNA wants a
# plain list of lists with a $data element and checks it with checkSets, so the
# structure is built with multiData() here and verified before the run.
#
# Set sizes are 4 to 8 libraries, which is small for consensus; minModuleSize is
# lowered to 20 and the outcome is reported together with checkSets output so the
# reader can see what the run actually had to work with.

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 16)

IN  <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_input"
OUT <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_out"
W <- 1000

M <- as.matrix(read.table(file.path(IN, sprintf("M_fib_w%d.tsv", W)), sep = "\t",
                          header = TRUE, row.names = 1, check.names = FALSE))
samples <- read.table(file.path(IN, sprintf("samples_w%d.tsv", W)), sep = "\t", header = TRUE)
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
  cs <- checkSets(multiExpr)
  cat(sprintf("после фильтра: генов %d, образцов: %s\n", cs$nGenes,
              paste(cs$nSamples, collapse = ", ")))
}

# WGCNA 1.74 has no multi-set soft-threshold picker, so the power is taken from
# the pooled fibroblast network of script 64, where pickSoftThreshold chose 14.
pw <- 14
cat(sprintf("беру softPower = %d, как в сети fib скрипта 64\n", pw))

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

me <- multiSetMEs(multiExpr, universalColors = cc)
for (i in seq_along(me)) {
  d <- data.frame(sample = rownames(multiExpr[[i]]$data), me[[i]]$data, check.names = FALSE)
  write.table(d, file.path(OUT, sprintf("consensus_ME_%s.tsv", names(multiExpr)[i])),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
cat("готово\n")
