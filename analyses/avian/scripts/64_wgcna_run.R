#!/usr/bin/env Rscript
# Standard WGCNA on the bird CAGE matrices built by script 62.
#
# Networks, all signed hybrid, all parameters at the package defaults recommended
# in Langfelder's tutorials:
#   fib      42 libraries, all seven species, no tissue contrast inside
#   epi      23 libraries, same
#   all_tc   65 libraries, every gene centred within tissue to drop the contrast
#   <species> 7-11 libraries each, below the package's own floor of 15, so these
#            are a sensitivity and are labelled as such
# plus blockwiseConsensusModules over the seven species as seven sets, inside
# fibroblasts, and modulePreservation of the fibroblast modules in the epiblast
# data and the other way round.
#
# Soft power: lowest value reaching scale-free R2 >= 0.8; if none does, the value
# from WGCNA's own recommendation table for the sample size (signed networks).

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 16)

IN  <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_input"
OUT <- "/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/wgcna_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
W <- 1000

read_m <- function(f) {
  m <- read.table(file.path(IN, f), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  as.matrix(m)
}
samples <- read.table(file.path(IN, sprintf("samples_w%d.tsv", W)), sep = "\t", header = TRUE)

pick_power <- function(datExpr, label) {
  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed hybrid",
                           verbose = 0)
  ok <- sft$fitIndices$SFT.R.sq >= 0.8
  n <- nrow(datExpr)
  fallback <- if (n < 20) 18 else if (n < 30) 16 else if (n < 40) 14 else 12
  p <- if (any(ok)) sft$fitIndices$Power[which(ok)[1]] else fallback
  write.table(cbind(label = label, sft$fitIndices), file.path(OUT, "soft_threshold.tsv"),
              sep = "\t", row.names = FALSE, append = file.exists(file.path(OUT, "soft_threshold.tsv")),
              col.names = !file.exists(file.path(OUT, "soft_threshold.tsv")), quote = FALSE)
  cat(sprintf("[%s] n=%d генов=%d   softPower=%d%s\n", label, n, ncol(datExpr), p,
              if (any(ok)) "" else "  (порог R2=0.8 не достигнут, значение из таблицы WGCNA)"))
  p
}

run_net <- function(m, label, minModuleSize = 30, center_by = NULL) {
  if (!is.null(center_by)) {
    for (lv in unique(center_by)) {
      j <- which(center_by == lv)
      m[, j] <- m[, j] - rowMeans(m[, j])
    }
  }
  datExpr <- t(m)
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    cat(sprintf("[%s] отброшено генов %d, образцов %d\n", label,
                sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
  }
  p <- pick_power(datExpr, label)
  net <- blockwiseModules(datExpr, power = p, networkType = "signed hybrid",
                          TOMType = "signed", minModuleSize = minModuleSize,
                          deepSplit = 2, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          maxBlockSize = 20000, verbose = 0)
  colors <- labels2colors(net$colors)
  me <- moduleEigengenes(datExpr, colors)$eigengenes
  kme <- signedKME(datExpr, me)
  gm <- data.frame(gene = colnames(datExpr), module = colors, label_num = net$colors,
                   stringsAsFactors = FALSE)
  best <- apply(kme, 1, function(r) max(r, na.rm = TRUE))
  gm$kME_own <- sapply(seq_len(nrow(gm)), function(i) {
    cn <- paste0("kME", gm$module[i])
    if (cn %in% colnames(kme)) kme[i, cn] else NA_real_
  })
  gm$kME_best <- best
  write.table(gm, file.path(OUT, sprintf("gene_modules_%s.tsv", label)), sep = "\t",
              row.names = FALSE, quote = FALSE)
  me2 <- data.frame(sample = rownames(datExpr), me, check.names = FALSE)
  write.table(me2, file.path(OUT, sprintf("module_eigengenes_%s.tsv", label)), sep = "\t",
              row.names = FALSE, quote = FALSE)
  sz <- as.data.frame(table(colors), stringsAsFactors = FALSE)
  names(sz) <- c("module", "n_genes")
  write.table(sz, file.path(OUT, sprintf("module_sizes_%s.tsv", label)), sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat(sprintf("[%s] модулей %d (без grey %d), grey %d генов\n", label,
              nrow(sz), sum(sz$module != "grey"),
              ifelse("grey" %in% sz$module, sz$n_genes[sz$module == "grey"], 0)))
  list(datExpr = datExpr, colors = colors, power = p)
}

cat("=== сети по тканям через все виды ===\n")
Mfib <- read_m(sprintf("M_fib_w%d.tsv", W))
Mepi <- read_m(sprintf("M_epi_w%d.tsv", W))
nf <- run_net(Mfib, "fib")
ne <- run_net(Mepi, "epi")

cat("\n=== обе ткани, центрирование внутри ткани ===\n")
Mall <- read_m(sprintf("M_all_w%d.tsv", W))
tis <- samples$tissue[match(colnames(Mall), samples$sample)]
na <- run_net(Mall, "all_tc", center_by = tis)

cat("\n=== по видам отдельно (чувствительность, образцов 7-11) ===\n")
for (sp in unique(samples$species)) {
  f <- file.path(IN, sprintf("M_%s_w%d.tsv", sp, W))
  if (!file.exists(f)) next
  m <- read_m(sprintf("M_%s_w%d.tsv", sp, W))
  keep <- apply(m, 1, function(r) sd(r) > 0)
  try(run_net(m[keep, , drop = FALSE], paste0("sp_", sp), minModuleSize = 20), silent = TRUE)
}

cat("\n=== consensus по семи видам как по семи sets, внутри фибробластов ===\n")
spp <- unique(samples$species)
multi <- vector("list", length(spp))
common <- rownames(Mfib)
for (i in seq_along(spp)) {
  cols <- grep(paste0("^", spp[i], "\\|fib\\|"), colnames(Mfib))
  multi[[i]] <- list(data = as.data.frame(t(Mfib[common, cols, drop = FALSE])))
}
names(multi) <- spp
multi <- structure(multi, class = "multiData")
pc <- min(sapply(multi, function(x) nrow(x$data)))
cat(sprintf("наборов %d, минимальный размер набора %d образцов\n", length(multi), pc))
cons <- try(blockwiseConsensusModules(multi, power = nf$power, networkType = "signed hybrid",
                                      TOMType = "signed", minModuleSize = 20, deepSplit = 2,
                                      mergeCutHeight = 0.25, numericLabels = TRUE,
                                      maxBlockSize = 20000, verbose = 0), silent = TRUE)
if (!inherits(cons, "try-error")) {
  cc <- labels2colors(cons$colors)
  write.table(data.frame(gene = common, module = cc), file.path(OUT, "consensus_modules_fib.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("consensus: модулей %d (без grey %d)\n", length(unique(cc)),
              length(setdiff(unique(cc), "grey"))))
} else {
  cat("consensus не сошёлся:", attr(cons, "condition")$message, "\n")
}

cat("\n=== сохранность модулей: фибробласты <-> эпибласт ===\n")
gcommon <- intersect(colnames(nf$datExpr), colnames(ne$datExpr))
cat(sprintf("общих генов %d\n", length(gcommon)))
mp_input <- list(fib = list(data = as.data.frame(nf$datExpr[, gcommon])),
                 epi = list(data = as.data.frame(ne$datExpr[, gcommon])))
col_fib <- nf$colors[match(gcommon, colnames(nf$datExpr))]
col_epi <- ne$colors[match(gcommon, colnames(ne$datExpr))]
mp <- try(modulePreservation(mp_input, list(fib = col_fib, epi = col_epi),
                             referenceNetworks = c(1, 2), nPermutations = 100,
                             randomSeed = 20260921, quickCor = 0, verbose = 0), silent = TRUE)
if (!inherits(mp, "try-error")) {
  for (ref in c(1, 2)) {
    test <- if (ref == 1) 2 else 1
    nm <- names(mp$preservation$Z[[ref]])[test]
    z <- mp$preservation$Z[[ref]][[test]]
    obs <- mp$preservation$observed[[ref]][[test]]
    out <- data.frame(module = rownames(z), moduleSize = obs$moduleSize,
                      Zsummary = z$Zsummary.pres, medianRank = obs$medianRank.pres)
    f <- file.path(OUT, sprintf("preservation_ref_%s.tsv", c("fib", "epi")[ref]))
    write.table(out, f, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("ссылка %s -> проверка %s: модулей %d, Zsummary>10 у %d, <2 у %d\n",
                c("fib", "epi")[ref], c("fib", "epi")[test], nrow(out),
                sum(out$Zsummary > 10, na.rm = TRUE), sum(out$Zsummary < 2, na.rm = TRUE)))
  }
} else {
  cat("modulePreservation не сошлась:", attr(mp, "condition")$message, "\n")
}
cat("\nготово, вывод в", OUT, "\n")
